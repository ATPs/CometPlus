// Copyright 2023 Jimmy Eng
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "NovelModeUtils.h"
#include <cctype>
#include <cstdio>
#include <fstream>

static bool StartsWithIgnoreCaseLocal(const string& sValue, const string& sPrefix)
{
   if (sValue.length() < sPrefix.length())
      return false;

   for (size_t i = 0; i < sPrefix.length(); ++i)
   {
      if (tolower((unsigned char)sValue[i]) != tolower((unsigned char)sPrefix[i]))
         return false;
   }

   return true;
}

static string SanitizeSourceLabelToken(const string& sInput)
{
   string sOut;
   sOut.reserve(sInput.length());

   for (size_t i = 0; i < sInput.length(); ++i)
   {
      unsigned char c = (unsigned char)sInput[i];
      if (isalnum(c))
      {
         sOut.push_back((char)c);
      }
      else
      {
         if (sOut.empty() || sOut[sOut.length() - 1] != '_')
            sOut.push_back('_');
      }
   }

   while (!sOut.empty() && sOut[sOut.length() - 1] == '_')
      sOut.erase(sOut.length() - 1);

   if (sOut.empty())
      sOut = "input";

   const size_t kMaxTokenLen = 48;
   if (sOut.length() > kMaxTokenLen)
      sOut.resize(kMaxTokenLen);

   return sOut;
}

bool BuildNovelMergedSourceLabel(size_t iOneBasedIndex,
                                 const string& sInputStem,
                                 string& sOutLabel,
                                 string& sErrorMsg)
{
   sOutLabel.clear();
   sErrorMsg.clear();

   if (iOneBasedIndex == 0)
   {
      sErrorMsg = " Error - invalid source label index for merged novel MGF.\n";
      return false;
   }

   string sToken = SanitizeSourceLabelToken(sInputStem);
   string sPrefix = "S" + std::to_string(iOneBasedIndex) + "_";
   sOutLabel = sPrefix + sToken;

   const size_t kMaxLabelLen = 64;
   if (sOutLabel.length() > kMaxLabelLen)
   {
      if (sPrefix.length() >= kMaxLabelLen)
         sOutLabel = sPrefix.substr(0, kMaxLabelLen);
      else
         sOutLabel = sPrefix + sToken.substr(0, kMaxLabelLen - sPrefix.length());
   }

   return true;
}

bool RewriteMgfTitleWithSourceLabel(const string& sTitleLine,
                                    const string& sSourceLabel,
                                    string& sOutTitleLine)
{
   sOutTitleLine = sTitleLine;
   if (!StartsWithIgnoreCaseLocal(sTitleLine, "TITLE="))
      return true;

   string sTitleValue = sTitleLine.substr(6);
   size_t iDot = sTitleValue.find('.');
   if (iDot != string::npos && iDot + 1 < sTitleValue.length())
      sOutTitleLine = "TITLE=" + sSourceLabel + sTitleValue.substr(iDot);
   else if (!sTitleValue.empty())
      sOutTitleLine = "TITLE=" + sSourceLabel + "." + sTitleValue;
   else
      sOutTitleLine = "TITLE=" + sSourceLabel;

   return true;
}

bool MergeFilteredMgfFilesWithSourceTag(
      const vector<std::pair<string, string>>& vFilteredMgfAndSourceLabel,
      string& sOutMergedMgfPath,
      string& sErrorMsg)
{
   sOutMergedMgfPath.clear();
   sErrorMsg.clear();

   if (vFilteredMgfAndSourceLabel.empty())
   {
      sErrorMsg = " Error - no filtered MGF files provided for merged novel search.\n";
      return false;
   }

   if (!CreateTempPath("cometplus_filtered_merged", ".mgf", sOutMergedMgfPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutMergedMgfPath.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create merged filtered MGF file \"" + sOutMergedMgfPath + "\".\n";
      sOutMergedMgfPath.clear();
      return false;
   }

   for (size_t i = 0; i < vFilteredMgfAndSourceLabel.size(); ++i)
   {
      const string& sInputMgf = vFilteredMgfAndSourceLabel.at(i).first;

      std::ifstream inFile(sInputMgf.c_str(), std::ios::in | std::ios::binary);
      if (!inFile.good())
      {
         outFile.close();
         remove(sOutMergedMgfPath.c_str());
         sErrorMsg = " Error - cannot read filtered MGF file \"" + sInputMgf + "\" for merge.\n";
         sOutMergedMgfPath.clear();
         return false;
      }

      string sLine;
      while (std::getline(inFile, sLine))
      {
         if (!sLine.empty() && sLine[sLine.length() - 1] == '\r')
            sLine.erase(sLine.length() - 1);

         outFile << sLine << "\n";
      }
   }

   outFile.flush();
   if (!outFile.good())
   {
      outFile.close();
      remove(sOutMergedMgfPath.c_str());
      sErrorMsg = " Error - failed while writing merged filtered MGF file \"" + sOutMergedMgfPath + "\".\n";
      sOutMergedMgfPath.clear();
      return false;
   }

   return true;
}

