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

#include "CometPlusInputFiles.h"
#include <cctype>
#include <fstream>

using namespace std;

static string TrimInputListLine(const string& sLine)
{
   size_t iBegin = 0;
   while (iBegin < sLine.length() && isspace((unsigned char)sLine[iBegin]))
      ++iBegin;

   size_t iEnd = sLine.length();
   while (iEnd > iBegin && isspace((unsigned char)sLine[iEnd - 1]))
      --iEnd;

   return sLine.substr(iBegin, iEnd - iBegin);
}

static bool ReadInputFilesList(const string& sPath,
                               vector<string>& vInputArgs,
                               string& sErrorMsg)
{
   vInputArgs.clear();
   sErrorMsg.clear();

   ifstream inFile(sPath.c_str(), ios::in);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot open --input_files list file \"" + sPath + "\".\n";
      return false;
   }

   string sLine;
   int iLineNumber = 0;
   while (getline(inFile, sLine))
   {
      ++iLineNumber;

      if (iLineNumber == 1
            && sLine.length() >= 3
            && (unsigned char)sLine[0] == 0xEF
            && (unsigned char)sLine[1] == 0xBB
            && (unsigned char)sLine[2] == 0xBF)
      {
         sLine = sLine.substr(3);
      }

      string sTrimmed = TrimInputListLine(sLine);
      if (sTrimmed.empty())
         continue;

      if (sTrimmed[0] == '#')
         continue;

      vInputArgs.push_back(sTrimmed);
   }

   if (inFile.bad())
   {
      sErrorMsg = " Error - failed while reading --input_files list file \"" + sPath + "\".\n";
      return false;
   }

   if (vInputArgs.empty())
   {
      sErrorMsg = " Error - --input_files list file \"" + sPath + "\" has no valid input paths.\n";
      return false;
   }

   return true;
}

bool ResolveInputArgsFromInputFilesOption(bool bInputFilesOptionSeen,
                                          const string& sInputFilesPath,
                                          vector<string>& vInputArgs,
                                          string& sErrorMsg)
{
   sErrorMsg.clear();

   if (!bInputFilesOptionSeen)
      return true;

   if (!vInputArgs.empty())
   {
      sErrorMsg = " Error - --input_files cannot be used together with positional input files.\n";
      return false;
   }

   vector<string> vInputArgsFromFile;
   if (!ReadInputFilesList(sInputFilesPath, vInputArgsFromFile, sErrorMsg))
      return false;

   vInputArgs.swap(vInputArgsFromFile);
   return true;
}
