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
#include "CometPlusMultiDB.h"
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

using namespace CometInterfaces;

extern comet_fileoffset_t clSizeCometFileOffset;
static string g_sCometPlusTempDir = ".";

bool IsIdxDatabasePath(const string& sPath)
{
   return (sPath.length() >= 4 && !STRCMP_IGNORE_CASE(sPath.c_str() + sPath.length() - 4, ".idx"));
}

static bool IsPathSeparator(char c)
{
   return c == '/' || c == '\\';
}

static bool IsAbsolutePathLocal(const string& sPath)
{
   if (sPath.empty())
      return false;

   if (IsPathSeparator(sPath[0]))
      return true;

#ifdef _WIN32
   if (sPath.length() >= 2
         && isalpha((unsigned char)sPath[0])
         && sPath[1] == ':')
      return true;
#endif

   return false;
}

static string JoinPathLocal(const string& sDir, const string& sName)
{
   if (sDir.empty() || sDir == ".")
      return sName;
   if (sName.empty())
      return sDir;
   if (IsAbsolutePathLocal(sName))
      return sName;
   if (IsPathSeparator(sDir[sDir.length() - 1]))
      return sDir + sName;
#ifdef _WIN32
   return sDir + "\\" + sName;
#else
   return sDir + "/" + sName;
#endif
}

void SetCometPlusTempDirectory(const string& sTempDir)
{
   if (sTempDir.empty())
      g_sCometPlusTempDir = ".";
   else
      g_sCometPlusTempDir = sTempDir;
}

bool EnsureDirectoryExistsRecursive(const string& sDirPath,
                                    string& sErrorMsg)
{
   sErrorMsg.clear();

   if (sDirPath.empty() || sDirPath == ".")
      return true;

   string sNorm = sDirPath;
   for (size_t i = 0; i < sNorm.length(); ++i)
   {
      if (sNorm[i] == '\\')
         sNorm[i] = '/';
   }

   string sCurrent;
   size_t iPos = 0;
#ifdef _WIN32
   if (sNorm.length() >= 2 && isalpha((unsigned char)sNorm[0]) && sNorm[1] == ':')
   {
      sCurrent = sNorm.substr(0, 2);
      iPos = 2;
   }
#endif
   if (iPos < sNorm.length() && sNorm[iPos] == '/')
   {
      sCurrent += "/";
      ++iPos;
   }

   while (iPos <= sNorm.length())
   {
      size_t iNext = sNorm.find('/', iPos);
      string sPart = (iNext == string::npos) ? sNorm.substr(iPos) : sNorm.substr(iPos, iNext - iPos);
      if (!sPart.empty() && sPart != ".")
      {
         if (!sCurrent.empty() && sCurrent[sCurrent.length() - 1] != '/')
            sCurrent += "/";
         sCurrent += sPart;

#ifdef _WIN32
         if (_mkdir(sCurrent.c_str()) != 0 && errno != EEXIST)
#else
         if (mkdir(sCurrent.c_str(), 0775) != 0 && errno != EEXIST)
#endif
         {
            struct stat st;
            if (stat(sCurrent.c_str(), &st) != 0 || !(st.st_mode & S_IFDIR))
            {
               sErrorMsg = " Error - cannot create directory \"" + sCurrent + "\".\n";
               return false;
            }
         }
      }

      if (iNext == string::npos)
         break;
      iPos = iNext + 1;
   }

   return true;
}


bool BuildMergedFasta(const vector<string>& vDatabases, string& sMergedDatabase, string& sErrorMsg)
{
   if (!CreateTempPath("cometplus_db", ".fasta", sMergedDatabase, sErrorMsg))
      return false;

   FILE* fpOut = fopen(sMergedDatabase.c_str(), "wb");
   if (fpOut == NULL)
   {
      sErrorMsg = " Error - cannot open temporary merged FASTA file \"" + sMergedDatabase + "\".\n";
      return false;
   }

   const size_t iBufSize = 1 << 20;
   vector<char> vBuf(iBufSize);

   for (size_t i = 0; i < vDatabases.size(); ++i)
   {
      FILE* fpIn = fopen(vDatabases.at(i).c_str(), "rb");
      if (fpIn == NULL)
      {
         fclose(fpOut);
         remove(sMergedDatabase.c_str());
         sMergedDatabase.clear();
         sErrorMsg = " Error - cannot read FASTA file \"" + vDatabases.at(i) + "\".\n";
         return false;
      }

      while (!feof(fpIn))
      {
         size_t iRead = fread(vBuf.data(), sizeof(char), vBuf.size(), fpIn);
         if (iRead > 0)
         {
            if (fwrite(vBuf.data(), sizeof(char), iRead, fpOut) != iRead)
            {
               fclose(fpIn);
               fclose(fpOut);
               remove(sMergedDatabase.c_str());
               sMergedDatabase.clear();
               sErrorMsg = " Error - failed while writing merged FASTA.\n";
               return false;
            }
         }
      }

      fclose(fpIn);
   }

   fclose(fpOut);
   return true;
}
static string TrimStringLocal(const string& sInput)
{
   size_t iStart = 0;
   while (iStart < sInput.length() && isspace((unsigned char)sInput[iStart]))
      ++iStart;

   size_t iEnd = sInput.length();
   while (iEnd > iStart && isspace((unsigned char)sInput[iEnd - 1]))
      --iEnd;

   return sInput.substr(iStart, iEnd - iStart);
}

static bool EndsWithIgnoreCase(const string& sValue, const string& sSuffix)
{
   if (sValue.length() < sSuffix.length())
      return false;

   size_t iOffset = sValue.length() - sSuffix.length();
   for (size_t i = 0; i < sSuffix.length(); ++i)
   {
      if (tolower((unsigned char)sValue[iOffset + i]) != tolower((unsigned char)sSuffix[i]))
         return false;
   }
   return true;
}

bool ConfigureDatabaseInputs(const vector<string>& vDatabases,
                             ICometSearchManager* pSearchMgr,
                             string& sMergedDatabasePath,
                             string& sErrorMsg)
{
   sMergedDatabasePath.clear();

   if (vDatabases.empty())
      return true;

   if (vDatabases.size() == 1)
   {
      pSearchMgr->SetParam("database_name", vDatabases.at(0), vDatabases.at(0));
      pSearchMgr->SetParam("database_name_list", "", string(""));
      pSearchMgr->SetParam("cometplus_multi_db_mode", "", string(""));
      return true;
   }

   bool bAllIdx = true;
   bool bAllFasta = true;
   for (size_t i = 0; i < vDatabases.size(); ++i)
   {
      bool bIdx = IsIdxDatabasePath(vDatabases.at(i));
      if (bIdx)
         bAllFasta = false;
      else
         bAllIdx = false;
   }

   if (!(bAllIdx || bAllFasta))
   {
      sErrorMsg = " Error - mixed database types are not allowed; provide either all FASTA files or all .idx files.\n";
      return false;
   }

   string sDbList;
   for (size_t i = 0; i < vDatabases.size(); ++i)
   {
      if (i > 0)
         sDbList += "\n";
      sDbList += vDatabases.at(i);
   }

   if (bAllIdx)
   {
      int iProbeType = 0;
      string sProbeError;
      for (size_t i = 0; i < vDatabases.size(); ++i)
      {
         int iObservedType = 0;
         if (!CometPlusProbeIdxType(vDatabases.at(i), iObservedType, sProbeError))
         {
            sErrorMsg = sProbeError;
            return false;
         }

         if (i == 0)
            iProbeType = iObservedType;
         else if (iObservedType != iProbeType)
         {
            sErrorMsg = " Error - mixed .idx types are not supported in one search invocation.\n";
            return false;
         }
      }

      pSearchMgr->SetParam("database_name", vDatabases.at(0), vDatabases.at(0));
      pSearchMgr->SetParam("database_name_list", sDbList, sDbList);
      pSearchMgr->SetParam("cometplus_multi_db_mode", "multi_idx", string("multi_idx"));
   }
   else
   {
      if (!BuildMergedFasta(vDatabases, sMergedDatabasePath, sErrorMsg))
         return false;

      pSearchMgr->SetParam("database_name", sMergedDatabasePath, sMergedDatabasePath);
      pSearchMgr->SetParam("database_name_list", sDbList, sDbList);
      pSearchMgr->SetParam("cometplus_multi_db_mode", "multi_fasta", string("multi_fasta"));
   }

   return true;
}

bool CreateTempPath(const string& sPrefix,
                    const string& sSuffix,
                    string& sOutPath,
                    string& sErrorMsg)
{
   string sTempRoot = g_sCometPlusTempDir.empty() ? "." : g_sCometPlusTempDir;
   if (!EnsureDirectoryExistsRecursive(sTempRoot, sErrorMsg))
      return false;

#ifdef _WIN32
   char szTmpPath[MAX_PATH];
   strncpy(szTmpPath, sTempRoot.c_str(), MAX_PATH - 1);
   szTmpPath[MAX_PATH - 1] = '\0';

   char szTmpFile[MAX_PATH];
   if (GetTempFileNameA(szTmpPath, "cpx", 0, szTmpFile) == 0)
   {
      sErrorMsg = " Error - cannot create temporary file name.\n";
      return false;
   }

   string sTmp = szTmpFile;
   if (!sSuffix.empty())
   {
      string sRenamed = sTmp + sSuffix;
      if (rename(sTmp.c_str(), sRenamed.c_str()) != 0)
      {
         sErrorMsg = " Error - cannot rename temporary file \"" + sTmp + "\".\n";
         return false;
      }
      sTmp = sRenamed;
   }

   sOutPath = sTmp;
   return true;
#else
   string sCleanPrefix = sPrefix;
   for (size_t i = 0; i < sCleanPrefix.length(); ++i)
   {
      if (!isalnum((unsigned char)sCleanPrefix[i]) && sCleanPrefix[i] != '_')
         sCleanPrefix[i] = '_';
   }
   if (sCleanPrefix.empty())
      sCleanPrefix = "cometplus";

   string sTemplate = JoinPathLocal(sTempRoot, sCleanPrefix + "_XXXXXX");
   vector<char> vTemplate(sTemplate.begin(), sTemplate.end());
   vTemplate.push_back('\0');

   int iFd = mkstemp(vTemplate.data());
   if (iFd == -1)
   {
      sErrorMsg = " Error - cannot create temporary file path.\n";
      return false;
   }
   close(iFd);

   string sTmpPath = vTemplate.data();
   if (!sSuffix.empty())
   {
      string sRenamed = sTmpPath + sSuffix;
      if (rename(sTmpPath.c_str(), sRenamed.c_str()) != 0)
      {
         remove(sTmpPath.c_str());
         sErrorMsg = " Error - cannot rename temporary file \"" + sTmpPath + "\".\n";
         return false;
      }
      sTmpPath = sRenamed;
   }

   sOutPath = sTmpPath;
   return true;
#endif
}

string NormalizePeptideToken(const string& sToken)
{
   string sOut;
   sOut.reserve(sToken.length());
   for (size_t i = 0; i < sToken.length(); ++i)
   {
      unsigned char c = (unsigned char)sToken[i];
      if (isalpha(c))
         sOut.push_back((char)toupper(c));
   }
   return sOut;
}

string NormalizePeptideForCompare(const string& sPeptide, bool bTreatSameIL)
{
   string sNorm = NormalizePeptideToken(sPeptide);
   if (bTreatSameIL)
   {
      for (size_t i = 0; i < sNorm.length(); ++i)
      {
         if (sNorm[i] == 'I')
            sNorm[i] = 'L';
      }
   }
   return sNorm;
}

string ComputeInputBaseName(const string& sInputPath)
{
   string sBase = sInputPath;
   size_t iSep = sBase.find_last_of("/\\");
   if (iSep != string::npos)
      sBase = sBase.substr(iSep + 1);

   if (EndsWithIgnoreCase(sBase, ".gz"))
      sBase = sBase.substr(0, sBase.length() - 3);

   size_t iDot = sBase.find_last_of('.');
   if (iDot != string::npos)
      sBase = sBase.substr(0, iDot);

   return sBase;
}

bool ParseScanIntegersFromString(const string& sInput,
                                 set<int>& setScans,
                                 string& sErrorMsg,
                                 const string& sSourceLabel)
{
   string sWork = sInput;
   for (size_t i = 0; i < sWork.length(); ++i)
   {
      if (sWork[i] == ',' || isspace((unsigned char)sWork[i]))
         sWork[i] = ' ';
   }

   std::istringstream iss(sWork);
   string sToken;
   while (iss >> sToken)
   {
      char* pEnd = NULL;
      long lVal = strtol(sToken.c_str(), &pEnd, 10);
      if (pEnd == NULL || *pEnd != '\0')
      {
         sErrorMsg = " Error - invalid scan token \"" + sToken + "\" in " + sSourceLabel + ".\n";
         return false;
      }
      if (lVal <= 0 || lVal > INT_MAX)
      {
         sErrorMsg = " Error - scan number out of range in " + sSourceLabel + ": \"" + sToken + "\".\n";
         return false;
      }
      setScans.insert((int)lVal);
   }

   return true;
}

bool ParseScanIntegersFromFile(const string& sPath,
                               set<int>& setScans,
                               string& sErrorMsg,
                               const string& sSourceLabel)
{
   std::ifstream inFile(sPath.c_str(), std::ios::in | std::ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read scan list file \"" + sPath + "\".\n";
      return false;
   }

   std::stringstream buffer;
   buffer << inFile.rdbuf();
   return ParseScanIntegersFromString(buffer.str(), setScans, sErrorMsg, sSourceLabel + " (" + sPath + ")");
}

bool ParseNovelPeptideFile(const string& sPath,
                           vector<string>& vPeptides,
                           string& sErrorMsg)
{
   std::ifstream inFile(sPath.c_str(), std::ios::in | std::ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read novel peptide file \"" + sPath + "\".\n";
      return false;
   }

   vector<string> vLines;
   string sLine;
   bool bFasta = false;
   while (std::getline(inFile, sLine))
   {
      vLines.push_back(sLine);
      string sTrimmed = TrimStringLocal(sLine);
      if (!sTrimmed.empty() && sTrimmed[0] == '>')
         bFasta = true;
   }

   unordered_set<string> setSeen;
   vPeptides.clear();

   if (bFasta)
   {
      string sCurrentSeq;
      for (size_t i = 0; i < vLines.size(); ++i)
      {
         string sTrimmed = TrimStringLocal(vLines.at(i));
         if (sTrimmed.empty())
            continue;

         if (sTrimmed[0] == '>')
         {
            if (!sCurrentSeq.empty())
            {
               string sNorm = NormalizePeptideToken(sCurrentSeq);
               if (!sNorm.empty() && setSeen.insert(sNorm).second)
                  vPeptides.push_back(sNorm);
               sCurrentSeq.clear();
            }
         }
         else
         {
            sCurrentSeq += sTrimmed;
         }
      }

      if (!sCurrentSeq.empty())
      {
         string sNorm = NormalizePeptideToken(sCurrentSeq);
         if (!sNorm.empty() && setSeen.insert(sNorm).second)
            vPeptides.push_back(sNorm);
      }
   }
   else
   {
      std::stringstream ss;
      for (size_t i = 0; i < vLines.size(); ++i)
      {
         if (i > 0)
            ss << "\n";
         ss << vLines.at(i);
      }

      string sRaw = ss.str();
      for (size_t i = 0; i < sRaw.length(); ++i)
      {
         if (sRaw[i] == ',' || isspace((unsigned char)sRaw[i]))
            sRaw[i] = ' ';
      }

      std::istringstream iss(sRaw);
      string sToken;
      while (iss >> sToken)
      {
         string sNorm = NormalizePeptideToken(sToken);
         if (!sNorm.empty() && setSeen.insert(sNorm).second)
            vPeptides.push_back(sNorm);
      }
   }

   return true;
}

static void SplitTabLine(const string& sLine, vector<string>& vFields)
{
   vFields.clear();
   size_t iStart = 0;
   while (iStart <= sLine.length())
   {
      size_t iPos = sLine.find('\t', iStart);
      if (iPos == string::npos)
      {
         vFields.push_back(sLine.substr(iStart));
         break;
      }

      vFields.push_back(sLine.substr(iStart, iPos - iStart));
      iStart = iPos + 1;
   }
}

bool ParseInternalNovelPeptideFile(const string& sPath,
                                   vector<NovelPeptideRecord>& vRecords,
                                   string& sErrorMsg)
{
   vRecords.clear();

   std::ifstream inFile(sPath.c_str(), std::ios::in | std::ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read internal novel peptide file \"" + sPath + "\".\n";
      return false;
   }

   string sHeader;
   if (!std::getline(inFile, sHeader))
   {
      sErrorMsg = " Error - internal novel peptide file is empty: \"" + sPath + "\".\n";
      return false;
   }

   vector<string> vHeaderFields;
   SplitTabLine(sHeader, vHeaderFields);

   int iPepCol = -1;
   int iPepIdCol = -1;
   int iProtIdCol = -1;
   for (size_t i = 0; i < vHeaderFields.size(); ++i)
   {
      string sCol = TrimStringLocal(vHeaderFields[i]);
      if (sCol == "peptide")
         iPepCol = (int)i;
      else if (sCol == "peptide_id")
         iPepIdCol = (int)i;
      else if (sCol == "protein_id")
         iProtIdCol = (int)i;
   }

   if (iPepCol < 0 || iPepIdCol < 0 || iProtIdCol < 0)
   {
      sErrorMsg = " Error - internal novel peptide file must contain header columns: peptide, peptide_id, protein_id.\n";
      return false;
   }

   unordered_map<string, size_t> mPepToIndex;
   string sLine;
   int iLineNumber = 1;
   while (std::getline(inFile, sLine))
   {
      ++iLineNumber;
      if (sLine.empty())
         continue;

      vector<string> vFields;
      SplitTabLine(sLine, vFields);
      int iNeeded = std::max(iPepCol, std::max(iPepIdCol, iProtIdCol));
      if ((int)vFields.size() <= iNeeded)
      {
         sErrorMsg = " Error - malformed internal novel peptide line " + std::to_string(iLineNumber) + ".\n";
         return false;
      }

      string sPeptide = NormalizePeptideToken(TrimStringLocal(vFields[iPepCol]));
      string sPeptideId = TrimStringLocal(vFields[iPepIdCol]);
      string sProteinField = TrimStringLocal(vFields[iProtIdCol]);

      if (sPeptide.empty() || sPeptideId.empty())
      {
         sErrorMsg = " Error - empty peptide or peptide_id at internal novel peptide line " + std::to_string(iLineNumber) + ".\n";
         return false;
      }

      vector<string> vProteinIds;
      unordered_set<string> setProteinSeen;
      std::stringstream ss(sProteinField);
      string sOneProtein;
      while (std::getline(ss, sOneProtein, ';'))
      {
         sOneProtein = TrimStringLocal(sOneProtein);
         if (!sOneProtein.empty() && setProteinSeen.insert(sOneProtein).second)
            vProteinIds.push_back(sOneProtein);
      }

      auto itExisting = mPepToIndex.find(sPeptide);
      if (itExisting == mPepToIndex.end())
      {
         NovelPeptideRecord rec;
         rec.sPeptide = sPeptide;
         rec.sPeptideId = sPeptideId;
         rec.vProteinIds = vProteinIds;
         mPepToIndex[sPeptide] = vRecords.size();
         vRecords.push_back(rec);
      }
      else
      {
         NovelPeptideRecord& rec = vRecords.at(itExisting->second);
         if (rec.sPeptideId != sPeptideId)
         {
            sErrorMsg = " Error - peptide \"" + sPeptide + "\" has conflicting peptide_id values in internal novel peptide file.\n";
            return false;
         }

         unordered_set<string> setMerge(rec.vProteinIds.begin(), rec.vProteinIds.end());
         for (size_t i = 0; i < vProteinIds.size(); ++i)
         {
            if (setMerge.insert(vProteinIds[i]).second)
               rec.vProteinIds.push_back(vProteinIds[i]);
         }
      }
   }

   if (vRecords.empty())
   {
      sErrorMsg = " Error - no peptide entries were parsed from internal novel peptide input.\n";
      return false;
   }

   return true;
}

bool FindNoCutEnzymeNumber(const string& sParamsFilePath,
                           int& iNoCutEnzyme,
                           string& sErrorMsg)
{
   std::ifstream inFile(sParamsFilePath.c_str());
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot open params file \"" + sParamsFilePath + "\" to resolve No_cut enzyme.\n";
      return false;
   }

   bool bInEnzymeSection = false;
   string sLine;
   while (std::getline(inFile, sLine))
   {
      string sTrimmed = TrimStringLocal(sLine);
      if (sTrimmed.empty())
         continue;

      if (!bInEnzymeSection)
      {
         if (sTrimmed == "[COMET_ENZYME_INFO]")
            bInEnzymeSection = true;
         continue;
      }

      if (sTrimmed[0] == '#')
         continue;

      int iNum = -1;
      char szName[128];
      if (sscanf(sTrimmed.c_str(), "%d. %127s", &iNum, szName) == 2)
      {
         string sName = szName;
         for (size_t i = 0; i < sName.length(); ++i)
            sName[i] = (char)tolower((unsigned char)sName[i]);

         if (sName == "no_cut" || sName == "nocut")
         {
            iNoCutEnzyme = iNum;
            return true;
         }
      }
   }

   sErrorMsg = " Error - could not resolve No_cut enzyme in params file \"" + sParamsFilePath + "\".\n";
   return false;
}

bool BuildTemporaryParamsFile(const string& sBaseParamsPath,
                              const map<string, string>& mOverrides,
                              string& sTmpParamsPath,
                              string& sErrorMsg)
{
   std::ifstream inFile(sBaseParamsPath.c_str(), std::ios::in | std::ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot open params file \"" + sBaseParamsPath + "\".\n";
      return false;
   }

   if (!CreateTempPath("cometplus_params", ".params", sTmpParamsPath, sErrorMsg))
      return false;

   std::ofstream outFile(sTmpParamsPath.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot write temporary params file \"" + sTmpParamsPath + "\".\n";
      return false;
   }

   map<string, string> mPending = mOverrides;
   bool bInserted = false;
   string sLine;
   while (std::getline(inFile, sLine))
   {
      string sRawLine = sLine;
      string sTrimmed = TrimStringLocal(sLine);

      if (!bInserted && sTrimmed == "[COMET_ENZYME_INFO]")
      {
         for (auto it = mPending.begin(); it != mPending.end(); ++it)
            outFile << it->first << " = " << it->second << "\n";
         bInserted = true;
      }

      bool bReplaced = false;
      size_t iEq = sLine.find('=');
      if (iEq != string::npos)
      {
         string sKey = TrimStringLocal(sLine.substr(0, iEq));
         if (!sKey.empty() && sKey[0] != '#')
         {
            auto itPending = mPending.find(sKey);
            if (itPending != mPending.end())
            {
               outFile << sKey << " = " << itPending->second << "\n";
               mPending.erase(itPending);
               bReplaced = true;
            }
         }
      }

      if (!bReplaced)
         outFile << sRawLine << "\n";
   }

   if (!bInserted)
   {
      for (auto it = mPending.begin(); it != mPending.end(); ++it)
         outFile << it->first << " = " << it->second << "\n";
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed while writing temporary params file \"" + sTmpParamsPath + "\".\n";
      return false;
   }

   return true;
}

static string EscapeShellArg(const string& sArg)
{
#ifdef _WIN32
   string sOut = "\"";
   for (size_t i = 0; i < sArg.length(); ++i)
   {
      if (sArg[i] == '"')
         sOut += "\\\"";
      else
         sOut += sArg[i];
   }
   sOut += "\"";
   return sOut;
#else
   string sOut = "'";
   for (size_t i = 0; i < sArg.length(); ++i)
   {
      if (sArg[i] == '\'')
         sOut += "'\"'\"'";
      else
         sOut += sArg[i];
   }
   sOut += "'";
   return sOut;
#endif
}

bool RunCometForIndexGeneration(const string& sExecutablePath,
                                const string& sParamsPath,
                                const string& sDatabasePath,
                                bool bFragmentIndex,
                                int iThreadOverride,
                                string& sErrorMsg)
{
   if (sExecutablePath.empty())
   {
      sErrorMsg = " Error - cannot resolve executable path for temporary index generation.\n";
      return false;
   }

   string sCmd = EscapeShellArg(sExecutablePath)
      + " --params " + EscapeShellArg(sParamsPath)
      + " --database " + EscapeShellArg(sDatabasePath)
      + " " + (bFragmentIndex ? "-i" : "-j");

   if (iThreadOverride > 0)
      sCmd += " --thread " + std::to_string(iThreadOverride);

   int iReturnCode = system(sCmd.c_str());
#ifdef _WIN32
   if (iReturnCode != 0)
   {
      sErrorMsg = " Error - temporary index generation failed for \"" + sDatabasePath + "\" with return code "
         + std::to_string(iReturnCode) + ".\n";
      return false;
   }
#else
   if (iReturnCode == -1)
   {
      sErrorMsg = " Error - failed to execute temporary index generation command.\n";
      return false;
   }
   if (!(WIFEXITED(iReturnCode) && WEXITSTATUS(iReturnCode) == 0))
   {
      sErrorMsg = " Error - temporary index generation failed for \"" + sDatabasePath + "\" with return code "
         + std::to_string(iReturnCode) + ".\n";
      return false;
   }
#endif

   return true;
}

bool ParsePeptideIdxEntries(const string& sIdxPath,
                            vector<PeptideMassEntry>& vEntries,
                            bool bLoadProteinIds,
                            string& sErrorMsg)
{
   vEntries.clear();
   FILE* fp = fopen(sIdxPath.c_str(), "rb");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot open peptide index file \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t lEndOfStruct = 0;
   comet_fileoffset_t lProteinsPos = 0;
   if (comet_fseek(fp, -clSizeCometFileOffset * 2, SEEK_END) != 0
         || fread(&lEndOfStruct, clSizeCometFileOffset, 1, fp) != 1
         || fread(&lProteinsPos, clSizeCometFileOffset, 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid peptide index footer in \"" + sIdxPath + "\".\n";
      return false;
   }

   vector<vector<string>> vProteinSets;
   if (bLoadProteinIds)
   {
      if (lProteinsPos <= 0 || comet_fseek(fp, lProteinsPos, SEEK_SET) != 0)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid peptide index protein section in \"" + sIdxPath + "\".\n";
         return false;
      }

      size_t tNumSets = 0;
      if (fread(&tNumSets, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - failed reading protein-set count from \"" + sIdxPath + "\".\n";
         return false;
      }

      vProteinSets.resize(tNumSets);
      for (size_t iSet = 0; iSet < tNumSets; ++iSet)
      {
         size_t tSetSize = 0;
         if (fread(&tSetSize, clSizeCometFileOffset, 1, fp) != 1)
         {
            fclose(fp);
            sErrorMsg = " Error - failed reading protein-set size from \"" + sIdxPath + "\".\n";
            return false;
         }

         vector<comet_fileoffset_t> vOffsets;
         vOffsets.resize(tSetSize);
         for (size_t i = 0; i < tSetSize; ++i)
         {
            if (fread(&vOffsets[i], clSizeCometFileOffset, 1, fp) != 1)
            {
               fclose(fp);
               sErrorMsg = " Error - failed reading protein offset from \"" + sIdxPath + "\".\n";
               return false;
            }
         }

         unordered_set<string> setSeen;
         for (size_t i = 0; i < vOffsets.size(); ++i)
         {
            if (comet_fseek(fp, vOffsets[i], SEEK_SET) != 0)
            {
               fclose(fp);
               sErrorMsg = " Error - failed seeking protein entry in \"" + sIdxPath + "\".\n";
               return false;
            }

            char szProteinName[512];
            if (fscanf(fp, "%500s", szProteinName) != 1)
            {
               fclose(fp);
               sErrorMsg = " Error - failed reading protein name from \"" + sIdxPath + "\".\n";
               return false;
            }
            szProteinName[500] = '\0';

            string sProteinId = szProteinName;
            if (!sProteinId.empty() && setSeen.insert(sProteinId).second)
               vProteinSets[iSet].push_back(sProteinId);
         }
      }
   }

   if (lEndOfStruct <= 0)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid peptide index structure offset in \"" + sIdxPath + "\".\n";
      return false;
   }

   if (comet_fseek(fp, lEndOfStruct, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - cannot seek peptide index structure in \"" + sIdxPath + "\".\n";
      return false;
   }

   int iMinMass = 0;
   int iMaxMass = 0;
   uint64_t tNumPeptides = 0;
   if (fread(&iMinMass, sizeof(int), 1, fp) != 1
         || fread(&iMaxMass, sizeof(int), 1, fp) != 1
         || fread(&tNumPeptides, sizeof(uint64_t), 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading peptide index summary for \"" + sIdxPath + "\".\n";
      return false;
   }
   (void)iMinMass;
   (void)tNumPeptides;

   if (iMaxMass <= 0)
   {
      fclose(fp);
      return true;
   }

   int iMaxMass10 = iMaxMass * 10;
   vector<comet_fileoffset_t> vIndex(iMaxMass10, -1);
   if (fread(vIndex.data(), clSizeCometFileOffset, iMaxMass10, fp) != (size_t)iMaxMass10)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading peptide index offsets from \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t lPeptideStart = -1;
   for (int i = 0; i < iMaxMass10; ++i)
   {
      if (vIndex[i] >= 0 && (lPeptideStart == -1 || vIndex[i] < lPeptideStart))
         lPeptideStart = vIndex[i];
   }
   if (lPeptideStart < 0)
   {
      fclose(fp);
      return true;
   }

   if (comet_fseek(fp, lPeptideStart, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - failed seeking peptide entries in \"" + sIdxPath + "\".\n";
      return false;
   }

   while (comet_ftell(fp) < lEndOfStruct)
   {
      int iLen = 0;
      if (fread(&iLen, sizeof(int), 1, fp) != 1)
         break;

      if (iLen <= 0 || iLen >= MAX_PEPTIDE_LEN)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid peptide length while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      string sPeptide;
      sPeptide.resize((size_t)iLen);
      if (fread(&sPeptide[0], sizeof(char), (size_t)iLen, fp) != (size_t)iLen)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide sequence while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      char cPrev = '-';
      char cNext = '-';
      if (fread(&cPrev, sizeof(char), 1, fp) != 1
            || fread(&cNext, sizeof(char), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide flanking residues while parsing \"" + sIdxPath + "\".\n";
         return false;
      }
      (void)cPrev;
      (void)cNext;

      unsigned char cNumMods = 0;
      if (fread(&cNumMods, sizeof(unsigned char), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide mod count while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      for (unsigned char iMod = 0; iMod < cNumMods; ++iMod)
      {
         unsigned char cPos = 0;
         char cWhichMod = 0;
         if (fread(&cPos, sizeof(unsigned char), 1, fp) != 1
               || fread(&cWhichMod, sizeof(char), 1, fp) != 1)
         {
            fclose(fp);
            sErrorMsg = " Error - truncated peptide mod sites while parsing \"" + sIdxPath + "\".\n";
            return false;
         }
         (void)cPos;
         (void)cWhichMod;
      }

      double dPepMass = 0.0;
      comet_fileoffset_t lProteinSetPos = 0;
      if (fread(&dPepMass, sizeof(double), 1, fp) != 1
            || fread(&lProteinSetPos, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide mass/protein-set entry while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      PeptideMassEntry entry;
      entry.sPeptide = sPeptide;
      entry.dMass = dPepMass;
      entry.bDecoy = false;
      if (bLoadProteinIds)
      {
         if (lProteinSetPos < 0 || (size_t)lProteinSetPos >= vProteinSets.size())
         {
            fclose(fp);
            sErrorMsg = " Error - invalid protein-set reference while parsing \"" + sIdxPath + "\".\n";
            return false;
         }

         entry.vProteinIds = vProteinSets[(size_t)lProteinSetPos];
      }
      vEntries.push_back(entry);
   }

   fclose(fp);
   return true;
}

bool ParseFragmentIdxEntries(const string& sIdxPath,
                             vector<PeptideMassEntry>& vEntries,
                             string& sErrorMsg)
{
   vEntries.clear();
   FILE* fp = fopen(sIdxPath.c_str(), "rb");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot open fragment index file \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t clPeptidesFilePos = 0;
   comet_fileoffset_t clProteinsFilePos = 0;
   comet_fileoffset_t clPermutationsFilePos = 0;
   if (comet_fseek(fp, -clSizeCometFileOffset * 3, SEEK_END) != 0
         || fread(&clPeptidesFilePos, clSizeCometFileOffset, 1, fp) != 1
         || fread(&clProteinsFilePos, clSizeCometFileOffset, 1, fp) != 1
         || fread(&clPermutationsFilePos, clSizeCometFileOffset, 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid fragment index footer in \"" + sIdxPath + "\".\n";
      return false;
   }
   (void)clProteinsFilePos;
   (void)clPermutationsFilePos;

   if (clPeptidesFilePos <= 0 || comet_fseek(fp, clPeptidesFilePos, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid fragment peptide section in \"" + sIdxPath + "\".\n";
      return false;
   }

   size_t tNumPeptides = 0;
   if (fread(&tNumPeptides, sizeof(size_t), 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading fragment peptide count from \"" + sIdxPath + "\".\n";
      return false;
   }

   for (size_t i = 0; i < tNumPeptides; ++i)
   {
      int iLen = 0;
      if (fread(&iLen, sizeof(int), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide entry while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      if (iLen <= 0 || iLen >= MAX_PEPTIDE_LEN)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid fragment peptide length while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      string sPeptide;
      sPeptide.resize((size_t)iLen);
      if (fread(&sPeptide[0], sizeof(char), (size_t)iLen, fp) != (size_t)iLen)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide sequence while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      char cPrev = '-';
      char cNext = '-';
      double dPepMass = 0.0;
      unsigned short siVarModProteinFilter = 0;
      comet_fileoffset_t lProteinSetPos = 0;

      if (fread(&cPrev, sizeof(char), 1, fp) != 1
            || fread(&cNext, sizeof(char), 1, fp) != 1
            || fread(&dPepMass, sizeof(double), 1, fp) != 1
            || fread(&siVarModProteinFilter, sizeof(unsigned short), 1, fp) != 1
            || fread(&lProteinSetPos, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide metadata while parsing \"" + sIdxPath + "\".\n";
         return false;
      }
      (void)cPrev;
      (void)cNext;
      (void)siVarModProteinFilter;
      (void)lProteinSetPos;

      PeptideMassEntry entry;
      entry.sPeptide = sPeptide;
      entry.dMass = dPepMass;
      entry.bDecoy = false;
      vEntries.push_back(entry);
   }

   fclose(fp);
   return true;
}

bool WritePeptidesToFasta(const vector<string>& vPeptides,
                          const string& sPrefix,
                          string& sOutPath,
                          string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".fasta", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   unordered_set<string> setSeen;
   int iWritten = 0;
   for (size_t i = 0; i < vPeptides.size(); ++i)
   {
      string sPep = NormalizePeptideToken(vPeptides.at(i));
      if (sPep.empty())
         continue;
      if (!setSeen.insert(sPep).second)
         continue;

      outFile << ">COMETPLUS_NOVEL_" << (iWritten + 1) << "\n";
      outFile << sPep << "\n";
      iWritten++;
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}

bool WriteNovelRecordsToFasta(const vector<NovelPeptideRecord>& vRecords,
                              const string& sPrefix,
                              string& sOutPath,
                              string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".fasta", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   unordered_set<string> setSeenPeptides;
   unordered_set<string> setSeenIds;
   int iAutoId = 1;

   for (size_t i = 0; i < vRecords.size(); ++i)
   {
      string sPep = NormalizePeptideToken(vRecords.at(i).sPeptide);
      if (sPep.empty())
         continue;
      if (!setSeenPeptides.insert(sPep).second)
         continue;

      string sPepId = TrimStringLocal(vRecords.at(i).sPeptideId);
      if (sPepId.empty())
      {
         sPepId = "COMETPLUS_NOVEL_" + std::to_string(iAutoId);
         iAutoId++;
      }

      if (!setSeenIds.insert(sPepId).second)
      {
         sErrorMsg = " Error - duplicate peptide_id encountered while building novel scoring FASTA: \"" + sPepId + "\".\n";
         return false;
      }

      outFile << ">" << sPepId << "\n";
      outFile << sPep << "\n";
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}

bool BuildNovelMassFilterContext(ICometSearchManager* pSearchMgr,
                                 NovelMassFilterContext& ctx,
                                 string& sErrorMsg)
{
   DoubleRange digestMassRange;
   if (!pSearchMgr->GetParamValue("digest_mass_range", digestMassRange))
   {
      sErrorMsg = " Error - could not read digest_mass_range parameter.\n";
      return false;
   }

   ctx.dPepMassLow = digestMassRange.dStart;
   ctx.dPepMassHigh = digestMassRange.dEnd;

   if (!pSearchMgr->GetParamValue("peptide_mass_tolerance_lower", ctx.dTolLower))
      ctx.dTolLower = -20.0;
   if (!pSearchMgr->GetParamValue("peptide_mass_tolerance_upper", ctx.dTolUpper))
      ctx.dTolUpper = 20.0;
   if (!pSearchMgr->GetParamValue("peptide_mass_units", ctx.iTolUnits))
      ctx.iTolUnits = 2;
   if (!pSearchMgr->GetParamValue("precursor_tolerance_type", ctx.iTolType))
      ctx.iTolType = 1;
   if (!pSearchMgr->GetParamValue("isotope_error", ctx.iIsotopeError))
      ctx.iIsotopeError = 0;
   if (!pSearchMgr->GetParamValue("mass_offsets", ctx.vMassOffsets))
      ctx.vMassOffsets.clear();

   IntRange precursorChargeRange;
   if (pSearchMgr->GetParamValue("precursor_charge", precursorChargeRange))
   {
      ctx.iStartCharge = precursorChargeRange.iStart;
      ctx.iEndCharge = precursorChargeRange.iEnd;
   }
   else
   {
      ctx.iStartCharge = 0;
      ctx.iEndCharge = 0;
   }

   if (!pSearchMgr->GetParamValue("override_charge", ctx.iOverrideCharge))
      ctx.iOverrideCharge = 0;
   if (!pSearchMgr->GetParamValue("min_precursor_charge", ctx.iMinPrecursorCharge))
      ctx.iMinPrecursorCharge = 1;
   if (!pSearchMgr->GetParamValue("max_precursor_charge", ctx.iMaxPrecursorCharge))
      ctx.iMaxPrecursorCharge = MAX_PRECURSOR_CHARGE;

   int iCorrectMass = 0;
   if (!pSearchMgr->GetParamValue("correct_mass", iCorrectMass))
      iCorrectMass = 0;
   ctx.bCorrectMass = (iCorrectMass != 0);

   if (!pSearchMgr->GetParamValue("ms_level", ctx.iMSLevel))
      ctx.iMSLevel = 2;

   return true;
}

static void BuildIsotopeShifts(const int iIsotopeError, vector<double>& vShifts)
{
   vShifts.clear();

   if (iIsotopeError == 7)
   {
      for (int x = -2; x <= 2; ++x)
         vShifts.push_back(x * 4.0070995);
      return;
   }

   int iMaxPositive = 0;
   if (iIsotopeError == 1)
      iMaxPositive = 1;
   else if (iIsotopeError == 2)
      iMaxPositive = 2;
   else if (iIsotopeError == 3 || iIsotopeError == 4 || iIsotopeError == 6)
      iMaxPositive = 3;
   else if (iIsotopeError == 5)
      iMaxPositive = 1;

   for (int x = 0; x <= iMaxPositive; ++x)
      vShifts.push_back((double)x * C13_DIFF);

   if (iIsotopeError == 4 || iIsotopeError == 5 || iIsotopeError == 6)
   {
      int iMaxNegative = (iIsotopeError == 6) ? 3 : 1;
      for (int x = 1; x <= iMaxNegative; ++x)
         vShifts.push_back(-(double)x * C13_DIFF);
   }

   if (vShifts.empty())
      vShifts.push_back(0.0);
}

static bool HasMassInRange(const vector<double>& vSortedMasses, double dLow, double dHigh)
{
   if (vSortedMasses.empty() || dLow > dHigh)
      return false;

   auto it = std::lower_bound(vSortedMasses.begin(), vSortedMasses.end(), dLow);
   if (it == vSortedMasses.end())
      return false;
   return (*it <= dHigh);
}

static bool MatchNovelMassForExperimental(const vector<double>& vSortedNovelMasses,
                                          double dExperimentalMass,
                                          int iCharge,
                                          const NovelMassFilterContext& ctx)
{
   if (vSortedNovelMasses.empty())
      return false;

   double dTolLow = ctx.dTolLower;
   double dTolHigh = ctx.dTolUpper;

   if (ctx.iTolUnits == 0) // amu
   {
      if (ctx.iTolType == 1)
      {
         dTolLow *= iCharge;
         dTolHigh *= iCharge;
      }
   }
   else if (ctx.iTolUnits == 1) // mmu
   {
      dTolLow *= 0.001;
      dTolHigh *= 0.001;
      if (ctx.iTolType == 1)
      {
         dTolLow *= iCharge;
         dTolHigh *= iCharge;
      }
   }
   else // ppm
   {
      double dMz = (dExperimentalMass + (iCharge - 1) * PROTON_MASS) / iCharge;
      double dBoundLower = dMz + dMz * ctx.dTolLower / 1.0E6;
      double dBoundUpper = dMz + dMz * ctx.dTolUpper / 1.0E6;
      double dProtonatedLower = dBoundLower * iCharge - iCharge * PROTON_MASS + PROTON_MASS;
      double dProtonatedUpper = dBoundUpper * iCharge - iCharge * PROTON_MASS + PROTON_MASS;
      dTolLow = dProtonatedLower - dExperimentalMass;
      dTolHigh = dProtonatedUpper - dExperimentalMass;
   }

   double dWindowLow = dExperimentalMass + dTolLow;
   double dWindowHigh = dExperimentalMass + dTolHigh;
   if (dWindowLow > dWindowHigh)
      std::swap(dWindowLow, dWindowHigh);

   vector<double> vShifts;
   BuildIsotopeShifts(ctx.iIsotopeError, vShifts);

   if (ctx.vMassOffsets.empty())
   {
      for (size_t i = 0; i < vShifts.size(); ++i)
      {
         if (HasMassInRange(vSortedNovelMasses, dWindowLow - vShifts[i], dWindowHigh - vShifts[i]))
            return true;
      }
   }
   else
   {
      for (size_t i = 0; i < ctx.vMassOffsets.size(); ++i)
      {
         for (size_t j = 0; j < vShifts.size(); ++j)
         {
            const double dDelta = ctx.vMassOffsets[i] + vShifts[j];
            if (HasMassInRange(vSortedNovelMasses, dWindowLow - dDelta, dWindowHigh - dDelta))
               return true;
         }
      }
   }

   return false;
}

bool FilterInputFileToTempMgf(const InputFileInfo& inputFile,
                              const set<int>& setExplicitScans,
                              bool bUseExplicitScans,
                              const vector<double>& vNovelMasses,
                              bool bUseNovelMassFilter,
                              const NovelMassFilterContext& ctx,
                              string& sOutTempMgfPath,
                              int& iNumScansKept,
                              string& sErrorMsg)
{
   iNumScansKept = 0;
   if (!CreateTempPath("cometplus_filtered_scans", ".mgf", sOutTempMgfPath, sErrorMsg))
      return false;

   FILE* fpOut = fopen(sOutTempMgfPath.c_str(), "w");
   if (fpOut == NULL)
   {
      sErrorMsg = " Error - cannot create temporary filtered spectrum file \"" + sOutTempMgfPath + "\".\n";
      return false;
   }

   vector<double> vSortedNovelMasses = vNovelMasses;
   std::sort(vSortedNovelMasses.begin(), vSortedNovelMasses.end());

   auto IsInRequestedRange = [&](int iScanNumber) -> bool
   {
      if (inputFile.iAnalysisType == AnalysisType_SpecificScan)
         return (iScanNumber == inputFile.iFirstScan);

      if (inputFile.iAnalysisType == AnalysisType_SpecificScanRange)
      {
         if (inputFile.iFirstScan > 0 && inputFile.iLastScan > 0)
            return (inputFile.iFirstScan <= iScanNumber && iScanNumber <= inputFile.iLastScan);
         if (inputFile.iFirstScan > 0)
            return (iScanNumber >= inputFile.iFirstScan);
         if (inputFile.iLastScan > 0)
            return (iScanNumber <= inputFile.iLastScan);
      }

      return true;
   };

   auto PassNovelMassPlausibility = [&](Spectrum& spec) -> bool
   {
      if (!bUseNovelMassFilter)
         return true;
      if (vSortedNovelMasses.empty())
         return false;

      int iSpectrumCharge = 0;
      for (int i = 0; i < spec.sizeMZ(); ++i)
      {
         double dMZ = 0.0;
         vector<pair<int, double>> vChargeStates;

         if (spec.sizeMZ() == spec.sizeZ())
         {
            iSpectrumCharge = spec.atZ(i).z;
         }
         else if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
         {
            iSpectrumCharge = spec.atZ(i).z;
         }
         else if (spec.sizeMZ() > spec.sizeZ())
         {
            iSpectrumCharge = 0;
         }
         else
         {
            iSpectrumCharge = 0;
         }

         dMZ = spec.getMonoMZ(i);

         if (ctx.bCorrectMass && spec.sizeMZ() == 1)
         {
            double dSelectionLower = spec.getSelWindowLower();
            double dSelectedMZ = spec.getMZ(i);
            if (dMZ > 0.1 && dSelectionLower > 0.1 && dMZ + 0.1 < dSelectionLower)
               dMZ = dSelectedMZ;
         }

         if (dMZ == 0.0)
            dMZ = spec.getMZ(i);

         if (dMZ == 0.0 && iSpectrumCharge != 0)
            dMZ = spec.atZ(i).mh / iSpectrumCharge;

         if (ctx.iStartCharge > 0 && ctx.iOverrideCharge > 0)
         {
            if (ctx.iOverrideCharge == 1)
            {
               for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
                  vChargeStates.push_back(make_pair(z, dMZ));
            }
            else if (ctx.iOverrideCharge == 2)
            {
               for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
               {
                  if (z == iSpectrumCharge)
                     vChargeStates.push_back(make_pair(z, dMZ));
               }
            }
            else if (ctx.iOverrideCharge == 3)
            {
               if (iSpectrumCharge > 0)
               {
                  vChargeStates.push_back(make_pair(iSpectrumCharge, dMZ));
               }
               else
               {
                  double dSumBelow = 0.0;
                  double dSumTotal = 0.0;
                  for (int ii = 0; ii < spec.size(); ++ii)
                  {
                     dSumTotal += spec.at(ii).intensity;
                     if (spec.at(ii).mz < spec.getMZ())
                        dSumBelow += spec.at(ii).intensity;
                  }

                  if (isEqual(dSumTotal, 0.0) || ((dSumBelow / dSumTotal) > 0.95))
                     vChargeStates.push_back(make_pair(1, dMZ));
                  else
                  {
                     for (int z = ctx.iStartCharge; z <= ctx.iEndCharge; ++z)
                        vChargeStates.push_back(make_pair(z, dMZ));
                  }
               }
            }
         }
         else
         {
            if (iSpectrumCharge > 0)
            {
               vChargeStates.push_back(make_pair(iSpectrumCharge, dMZ));

               if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
               {
                  for (int ii = 1; ii < spec.sizeZ(); ++ii)
                  {
                     vChargeStates.push_back(make_pair(spec.atZ(ii).z,
                        (spec.atZ(ii).mh + PROTON_MASS * (spec.atZ(ii).z - 1)) / spec.atZ(ii).z));
                  }
               }
            }
            else
            {
               double dSumBelow = 0.0;
               double dSumTotal = 0.0;
               for (int ii = 0; ii < spec.size(); ++ii)
               {
                  dSumTotal += spec.at(ii).intensity;
                  if (spec.at(ii).mz < spec.getMZ())
                     dSumBelow += spec.at(ii).intensity;
               }

               if (isEqual(dSumTotal, 0.0) || ((dSumBelow / dSumTotal) > 0.95))
               {
                  vChargeStates.push_back(make_pair(1, dMZ));
               }
               else
               {
                  vChargeStates.push_back(make_pair(2, dMZ));
                  vChargeStates.push_back(make_pair(3, dMZ));
               }
            }
         }

         for (size_t iz = 0; iz < vChargeStates.size(); ++iz)
         {
            int iPrecursorCharge = vChargeStates[iz].first;
            double dExperimentalMass = vChargeStates[iz].second * iPrecursorCharge
               - PROTON_MASS * (iPrecursorCharge - 1);

            if ((isEqual(ctx.dPepMassLow, 0.0)
                     || (dExperimentalMass >= ctx.dPepMassLow && dExperimentalMass <= ctx.dPepMassHigh))
                  && iPrecursorCharge <= ctx.iMaxPrecursorCharge
                  && iPrecursorCharge >= ctx.iMinPrecursorCharge)
            {
               if (MatchNovelMassForExperimental(vSortedNovelMasses, dExperimentalMass, iPrecursorCharge, ctx))
                  return true;
            }
         }
      }

      return false;
   };

   MSReader mstReader;
   vector<MSSpectrumType> msLevelFilter;
   if (ctx.iMSLevel == 1)
      msLevelFilter.push_back(MS1);
   else if (ctx.iMSLevel == 3)
      msLevelFilter.push_back(MS3);
   else
      msLevelFilter.push_back(MS2);
   mstReader.setFilter(msLevelFilter);

   Spectrum spec;
   bool bRead = mstReader.readFile(inputFile.szFileName, spec, 0);
   while (bRead)
   {
      int iScanNumber = spec.getScanNumber();

      if (iScanNumber > 0 && IsInRequestedRange(iScanNumber))
      {
         bool bPassExplicit = true;
         if (bUseExplicitScans)
            bPassExplicit = (setExplicitScans.find(iScanNumber) != setExplicitScans.end());

         if (bPassExplicit && PassNovelMassPlausibility(spec))
         {
            mstReader.appendMGFFileHandle(fpOut, spec);
            iNumScansKept++;
         }
      }

      bRead = mstReader.readFile(NULL, spec);
   }

   fclose(fpOut);
   return true;
}
