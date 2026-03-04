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

#include "CometPlusRuntimeUtils.h"
#include <cctype>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <map>
#include <unordered_set>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

using namespace std;

static bool IsPathSeparatorLocal(char c)
{
   return c == '/' || c == '\\';
}

bool ContainsPathSeparator(const string& sPath)
{
   for (size_t i = 0; i < sPath.length(); ++i)
   {
      if (IsPathSeparatorLocal(sPath[i]))
         return true;
   }
   return false;
}

bool IsAbsolutePathLocal(const string& sPath)
{
   if (sPath.empty())
      return false;

   if (IsPathSeparatorLocal(sPath[0]))
      return true;

#ifdef _WIN32
   if (sPath.length() >= 2 && isalpha((unsigned char)sPath[0]) && sPath[1] == ':')
      return true;
#endif

   return false;
}

string JoinPathLocal(const string& sDir, const string& sName)
{
   if (sDir.empty() || sDir == ".")
      return sName;
   if (sName.empty())
      return sDir;
   if (IsAbsolutePathLocal(sName))
      return sName;
   if (IsPathSeparatorLocal(sDir[sDir.length() - 1]))
      return sDir + sName;
#ifdef _WIN32
   return sDir + "\\" + sName;
#else
   return sDir + "/" + sName;
#endif
}

string GetDirNameLocal(const string& sPath)
{
   size_t iPos = sPath.find_last_of("/\\");
   if (iPos == string::npos)
      return "";
   if (iPos == 0)
      return sPath.substr(0, 1);
   return sPath.substr(0, iPos);
}

bool FileExistsLocal(const string& sPath)
{
   struct stat st;
   return stat(sPath.c_str(), &st) == 0;
}

string NormalizePathKeyLocal(const string& sPath)
{
   string sNorm = sPath;
   for (size_t i = 0; i < sNorm.length(); ++i)
   {
      if (sNorm[i] == '\\')
         sNorm[i] = '/';
   }

   while (sNorm.length() > 1 && sNorm[sNorm.length() - 1] == '/')
      sNorm.erase(sNorm.length() - 1);

   return sNorm;
}

bool EndsWithIgnoreCaseLocal(const string& sValue, const string& sSuffix)
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

bool IsMzMLbPath(const string& sPath)
{
   return EndsWithIgnoreCaseLocal(sPath, ".mzMLb");
}

bool HasMzMLbInputs(const vector<InputFileInfo*>& vInputs)
{
   for (size_t i = 0; i < vInputs.size(); ++i)
   {
      if (vInputs.at(i) != NULL && IsMzMLbPath(vInputs.at(i)->szFileName))
         return true;
   }

   return false;
}

bool IsExecutableFileLocal(const string& sPath)
{
   if (sPath.empty())
      return false;

   struct stat st;
   if (stat(sPath.c_str(), &st) != 0)
      return false;

#ifdef _WIN32
   return true;
#else
   if (!S_ISREG(st.st_mode))
      return false;
   return access(sPath.c_str(), X_OK) == 0;
#endif
}

bool FindExecutableOnPath(const string& sExeName, string& sOutPath)
{
   sOutPath.clear();
   if (sExeName.empty())
      return false;

   if (ContainsPathSeparator(sExeName) || IsAbsolutePathLocal(sExeName))
   {
      if (IsExecutableFileLocal(sExeName))
      {
         sOutPath = sExeName;
         return true;
      }
      return false;
   }

   const char* pszPath = getenv("PATH");
   if (pszPath == NULL || *pszPath == '\0')
      return false;

   const string sPathValue = pszPath;
#ifdef _WIN32
   const char cSep = ';';
#else
   const char cSep = ':';
#endif

   size_t iStart = 0;
   while (iStart <= sPathValue.length())
   {
      size_t iSep = sPathValue.find(cSep, iStart);
      string sDir;
      if (iSep == string::npos)
      {
         sDir = sPathValue.substr(iStart);
         iStart = sPathValue.length() + 1;
      }
      else
      {
         sDir = sPathValue.substr(iStart, iSep - iStart);
         iStart = iSep + 1;
      }

      if (sDir.empty())
         sDir = ".";

      string sCandidate = JoinPathLocal(sDir, sExeName);
      if (IsExecutableFileLocal(sCandidate))
      {
         sOutPath = sCandidate;
         return true;
      }
   }

   return false;
}

bool ResolvePrefilterWorkerExecutablePath(const string& sMainExePath,
                                          string& sOutPath,
                                          string& sErrorMsg)
{
   sOutPath.clear();
   sErrorMsg.clear();

   const char* pszOverride = getenv("COMETPLUS_PREFILTER_WORKER");
   if (pszOverride != NULL && *pszOverride != '\0')
   {
      string sOverride = pszOverride;
      string sFoundOnPath;

      if (ContainsPathSeparator(sOverride) || IsAbsolutePathLocal(sOverride))
      {
         if (IsExecutableFileLocal(sOverride))
         {
            sOutPath = sOverride;
            return true;
         }

         sErrorMsg = " Error - COMETPLUS_PREFILTER_WORKER points to a non-executable path: \""
            + sOverride + "\".\n";
         return false;
      }

      if (FindExecutableOnPath(sOverride, sFoundOnPath))
      {
         sOutPath = sFoundOnPath;
         return true;
      }

      if (IsExecutableFileLocal(sOverride))
      {
         sOutPath = sOverride;
         return true;
      }

      sErrorMsg = " Error - COMETPLUS_PREFILTER_WORKER=\"" + sOverride
         + "\" was not found in PATH and is not executable in current directory.\n";
      return false;
   }

   string sSiblingPath;
   string sExeDir = GetDirNameLocal(sMainExePath);
   if (!sExeDir.empty())
   {
      sSiblingPath = JoinPathLocal(sExeDir, "cometplus_prefilter_worker");
      if (IsExecutableFileLocal(sSiblingPath))
      {
         sOutPath = sSiblingPath;
         return true;
      }
   }

   string sFoundOnPath;
   if (FindExecutableOnPath("cometplus_prefilter_worker", sFoundOnPath))
   {
      sOutPath = sFoundOnPath;
      return true;
   }

   if (IsExecutableFileLocal("cometplus_prefilter_worker"))
   {
      sOutPath = "cometplus_prefilter_worker";
      return true;
   }

   sErrorMsg = string(" Error - mzMLb parallel prefilter worker executable not found.")
      + " Checked sibling path \"" + (sSiblingPath.empty() ? "<unresolved>" : sSiblingPath) + "\""
      + " and searched PATH for \"cometplus_prefilter_worker\"."
      + " You can override via COMETPLUS_PREFILTER_WORKER=/abs/path/to/cometplus_prefilter_worker.\n";
   return false;
}

static string EscapeShellArgLocal(const string& sArg)
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

bool WriteIntegerSetToTempFile(const set<int>& setValues,
                               const string& sPrefix,
                               string& sOutPath,
                               string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".txt", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary list file \"" + sOutPath + "\".\n";
      return false;
   }

   for (auto it = setValues.begin(); it != setValues.end(); ++it)
      outFile << *it << "\n";

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary list file \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}

bool WriteDoubleVectorToTempFile(const vector<double>& vValues,
                                 const string& sPrefix,
                                 string& sOutPath,
                                 string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".txt", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary list file \"" + sOutPath + "\".\n";
      return false;
   }

   char szBuf[64];
   for (size_t i = 0; i < vValues.size(); ++i)
   {
      snprintf(szBuf, sizeof(szBuf), "%.17g", vValues.at(i));
      outFile << szBuf << "\n";
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary list file \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}

bool WritePrefilterWorkerJobFile(const InputFileInfo& inputFile,
                                 bool bUseExplicitScans,
                                 const string& sExplicitScansFile,
                                 bool bUseNovelMassFilter,
                                 const string& sNovelMassesFile,
                                 const string& sInternalNovelPeptideFile,
                                 const NovelMassFilterContext& ctx,
                                 const string& sMassOffsetsFile,
                                 const string& sTempDir,
                                 const string& sResultFile,
                                 string& sOutJobFile,
                                 string& sErrorMsg)
{
   if (!CreateTempPath("cometplus_prefilter_job", ".tsv", sOutJobFile, sErrorMsg))
      return false;

   std::ofstream outFile(sOutJobFile.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create prefilter worker job file \"" + sOutJobFile + "\".\n";
      return false;
   }

   auto WriteKV = [&](const string& sKey, const string& sValue)
   {
      outFile << sKey << "\t" << sValue << "\n";
   };

   auto BoolTo01 = [&](bool bValue) -> string
   {
      return bValue ? "1" : "0";
   };

   char szBuf[64];
   WriteKV("input_file", inputFile.szFileName);
   WriteKV("analysis_type", std::to_string(inputFile.iAnalysisType));
   WriteKV("first_scan", std::to_string(inputFile.iFirstScan));
   WriteKV("last_scan", std::to_string(inputFile.iLastScan));
   WriteKV("use_explicit_scans", BoolTo01(bUseExplicitScans));
   WriteKV("explicit_scans_file", sExplicitScansFile);
   WriteKV("use_novel_mass_filter", BoolTo01(bUseNovelMassFilter));
   WriteKV("novel_masses_file", sNovelMassesFile);
   WriteKV("internal_novel_peptide_file", sInternalNovelPeptideFile);
   WriteKV("mass_offsets_file", sMassOffsetsFile);
   WriteKV("temp_dir", sTempDir);
   WriteKV("result_file", sResultFile);

   snprintf(szBuf, sizeof(szBuf), "%.17g", ctx.dPepMassLow);
   WriteKV("ctx_dPepMassLow", szBuf);
   snprintf(szBuf, sizeof(szBuf), "%.17g", ctx.dPepMassHigh);
   WriteKV("ctx_dPepMassHigh", szBuf);
   snprintf(szBuf, sizeof(szBuf), "%.17g", ctx.dTolLower);
   WriteKV("ctx_dTolLower", szBuf);
   snprintf(szBuf, sizeof(szBuf), "%.17g", ctx.dTolUpper);
   WriteKV("ctx_dTolUpper", szBuf);
   WriteKV("ctx_iTolUnits", std::to_string(ctx.iTolUnits));
   WriteKV("ctx_iTolType", std::to_string(ctx.iTolType));
   WriteKV("ctx_iIsotopeError", std::to_string(ctx.iIsotopeError));
   WriteKV("ctx_iStartCharge", std::to_string(ctx.iStartCharge));
   WriteKV("ctx_iEndCharge", std::to_string(ctx.iEndCharge));
   WriteKV("ctx_iOverrideCharge", std::to_string(ctx.iOverrideCharge));
   WriteKV("ctx_iMinPrecursorCharge", std::to_string(ctx.iMinPrecursorCharge));
   WriteKV("ctx_iMaxPrecursorCharge", std::to_string(ctx.iMaxPrecursorCharge));
   WriteKV("ctx_bCorrectMass", BoolTo01(ctx.bCorrectMass));
   WriteKV("ctx_iMSLevel", std::to_string(ctx.iMSLevel));

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing prefilter worker job file \"" + sOutJobFile + "\".\n";
      return false;
   }

   return true;
}

static bool ParseIntStrictLocal(const string& sValue, int& iOut)
{
   char* pEnd = NULL;
   long lValue = strtol(sValue.c_str(), &pEnd, 10);
   if (pEnd == NULL || *pEnd != '\0')
      return false;
   if (lValue < INT_MIN || lValue > INT_MAX)
      return false;
   iOut = (int)lValue;
   return true;
}

bool ParsePrefilterWorkerResultFile(const string& sResultFile,
                                    PrefilterWorkerOutput& output,
                                    string& sErrorMsg)
{
   output = PrefilterWorkerOutput();

   std::ifstream inFile(sResultFile.c_str(), std::ios::in | std::ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read prefilter worker result file \"" + sResultFile + "\".\n";
      return false;
   }

   std::map<string, string> mValues;
   string sLine;
   while (std::getline(inFile, sLine))
   {
      if (!sLine.empty() && sLine[sLine.length() - 1] == '\r')
         sLine.erase(sLine.length() - 1);
      if (sLine.empty())
         continue;

      size_t iTab = sLine.find('\t');
      if (iTab == string::npos)
         continue;

      mValues[sLine.substr(0, iTab)] = sLine.substr(iTab + 1);
   }

   auto itSuccess = mValues.find("success");
   if (itSuccess == mValues.end())
   {
      sErrorMsg = " Error - malformed prefilter worker result file \"" + sResultFile + "\" (missing success).\n";
      return false;
   }

   int iSuccess = 0;
   if (!ParseIntStrictLocal(itSuccess->second, iSuccess) || (iSuccess != 0 && iSuccess != 1))
   {
      sErrorMsg = " Error - malformed prefilter worker result file \"" + sResultFile + "\" (invalid success).\n";
      return false;
   }
   output.bSuccess = (iSuccess == 1);

   auto itPath = mValues.find("temp_mgf_path");
   if (itPath != mValues.end())
      output.sTempMgfPath = itPath->second;

   auto itCount = mValues.find("num_scans_kept");
   if (itCount != mValues.end())
   {
      int iCount = 0;
      if (ParseIntStrictLocal(itCount->second, iCount))
         output.iNumScansKept = iCount;
   }

   auto itError = mValues.find("error");
   if (itError != mValues.end())
      output.sErrorMsg = itError->second;

   return true;
}

bool RunPrefilterWorkerJob(const string& sWorkerExePath,
                           const string& sJobFile,
                           string& sErrorMsg)
{
   string sCmd = EscapeShellArgLocal(sWorkerExePath)
      + " --job " + EscapeShellArgLocal(sJobFile);

   int iReturnCode = system(sCmd.c_str());
#ifdef _WIN32
   if (iReturnCode != 0)
   {
      sErrorMsg = " Error - prefilter worker command failed with return code "
         + std::to_string(iReturnCode) + ".\n";
      return false;
   }
#else
   if (iReturnCode == -1)
   {
      sErrorMsg = " Error - failed to execute prefilter worker command.\n";
      return false;
   }
   if (!(WIFEXITED(iReturnCode) && WEXITSTATUS(iReturnCode) == 0))
   {
      sErrorMsg = " Error - prefilter worker command failed with return code "
         + std::to_string(iReturnCode) + ".\n";
      return false;
   }
#endif

   return true;
}

string ResolveInternalOutputPath(const string& sPath, const string& sOutputFolder)
{
   if (sPath.empty())
      return sPath;

   if (IsAbsolutePathLocal(sPath) || ContainsPathSeparator(sPath))
      return sPath;

   return JoinPathLocal(sOutputFolder, sPath);
}

string GetLocalTimestampString()
{
   time_t tNow;
   time(&tNow);
   char szTime[64];
   szTime[0] = '\0';
   strftime(szTime, sizeof(szTime), "%Y-%m-%d %H:%M:%S", localtime(&tNow));
   return szTime;
}

const char* BoolToOnOffString(bool bValue)
{
   return bValue ? "on" : "off";
}

void LogRetainedTempArtifacts(const vector<string>& vTempArtifacts,
                              const string& sMergedDatabasePath)
{
   vector<string> vRetained;
   unordered_set<string> setSeen;

   auto AddIfUnique = [&](const string& sPath)
   {
      if (sPath.empty())
         return;
      if (!setSeen.insert(sPath).second)
         return;
      vRetained.push_back(sPath);
   };

   AddIfUnique(sMergedDatabasePath);
   for (size_t i = 0; i < vTempArtifacts.size(); ++i)
      AddIfUnique(vTempArtifacts.at(i));

   char szBuf[1024];
   if (vRetained.empty())
   {
      snprintf(szBuf,
               sizeof(szBuf),
               " [%s] --keep-tmp enabled: no temporary artifacts were recorded.\n",
               GetLocalTimestampString().c_str());
      logout(szBuf);
      return;
   }

   snprintf(szBuf,
            sizeof(szBuf),
            " [%s] --keep-tmp enabled: retained %zu temporary artifacts for debugging.\n",
            GetLocalTimestampString().c_str(),
            vRetained.size());
   logout(szBuf);

   for (size_t i = 0; i < vRetained.size(); ++i)
   {
      snprintf(szBuf, sizeof(szBuf), "   - %s\n", vRetained.at(i).c_str());
      logout(szBuf);
   }
}

void LogStageTiming(const string& sStage,
                    const std::chrono::steady_clock::time_point& tStageStart,
                    const std::chrono::steady_clock::time_point& tProgramStart)
{
   const auto tNow = std::chrono::steady_clock::now();
   double dStageSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(tNow - tStageStart).count() / 1000.0;
   double dTotalSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(tNow - tProgramStart).count() / 1000.0;

   char szBuf[1024];
   snprintf(szBuf,
            sizeof(szBuf),
            " [%s] %s done (%.3f sec; total %.3f sec)\n",
            GetLocalTimestampString().c_str(),
            sStage.c_str(),
            dStageSeconds,
            dTotalSeconds);
   logout(szBuf);
}

bool PathHasSeparatorForName(const string& sName)
{
   return ContainsPathSeparator(sName);
}

void AppendPlannedOutputFilesForInput(const InputFileInfo& inputInfo,
                                      const string& sBaseName,
                                      const string& sOutputSuffix,
                                      const string& sTxtExt,
                                      bool bOutputSqtFile,
                                      bool bOutputTxtFile,
                                      bool bOutputPepXmlFile,
                                      int iOutputMzidFile,
                                      bool bOutputPercolatorFile,
                                      int iDecoySearch,
                                      vector<string>& vPlanned)
{
   auto AddPairSuffix = [&](const string& sDefaultExt,
                            const string& sTargetExt,
                            const string& sDecoyExt)
   {
#ifndef CRUX
      (void)sTargetExt;
#endif
      bool bEntire = (inputInfo.iAnalysisType == AnalysisType_EntireFile);
      string sRangeToken;
      if (!bEntire)
      {
         sRangeToken = "." + std::to_string(inputInfo.iFirstScan) + "-" + std::to_string(inputInfo.iLastScan);
      }

      string sMain;
      if (bEntire)
         sMain = sBaseName + sOutputSuffix + sDefaultExt;
      else
         sMain = sBaseName + sOutputSuffix + sRangeToken + sDefaultExt;

#ifdef CRUX
      if (iDecoySearch == 2)
      {
         if (bEntire)
            sMain = sBaseName + sOutputSuffix + sTargetExt;
         else
            sMain = sBaseName + sOutputSuffix + sRangeToken + sTargetExt;
      }
#endif

      vPlanned.push_back(sMain);

      if (iDecoySearch == 2)
      {
         string sDecoy;
         if (bEntire)
            sDecoy = sBaseName + sOutputSuffix + sDecoyExt;
         else
            sDecoy = sBaseName + sOutputSuffix + sRangeToken + sDecoyExt;
         vPlanned.push_back(sDecoy);
      }
   };

   if (bOutputSqtFile)
      AddPairSuffix(".sqt", ".target.sqt", ".decoy.sqt");

   if (bOutputTxtFile)
      AddPairSuffix("." + sTxtExt, ".target." + sTxtExt, ".decoy." + sTxtExt);

   if (bOutputPepXmlFile)
      AddPairSuffix(".pep.xml", ".target.pep.xml", ".decoy.pep.xml");

   if (iOutputMzidFile != 0)
      AddPairSuffix(".mzid", ".target.mzid", ".decoy.mzid");

   if (bOutputPercolatorFile)
   {
      bool bEntire = (inputInfo.iAnalysisType == AnalysisType_EntireFile);
      if (bEntire)
         vPlanned.push_back(sBaseName + sOutputSuffix + ".pin");
      else
      {
         vPlanned.push_back(sBaseName + sOutputSuffix + "."
               + std::to_string(inputInfo.iFirstScan) + "-" + std::to_string(inputInfo.iLastScan) + ".pin");
      }
   }
}
