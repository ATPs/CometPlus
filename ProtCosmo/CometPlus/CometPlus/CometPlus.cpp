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


#include "Common.h"
#include "CometData.h"
#include "CometDataInternal.h"
#include "CometInterfaces.h"
#include "CometPlusMultiDB.h"
#include "CometPlusParams.h"
#include "NovelModeUtils.h"
#include "githubsha.h"
#include <algorithm>
#include <chrono>
#include <cerrno>
#include <climits>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <functional>
#include <atomic>
#include <mutex>
#include <thread>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

using namespace CometInterfaces;

static string g_sCometPlusExecutablePath;
static bool g_bCometPlusKeepTempFiles = false;

void Usage(char *pszCmd,
           bool bFullHelp = false,
           const char *pszParamsFile = NULL);
void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    ICometSearchManager *pSearchMgr,
                    string &sMergedDatabasePath,
                    vector<string> &vTempArtifacts);
void PrintParams(int iPrintParams);
bool ValidateInputFile(char *pszInputFileName);
static void LogRetainedTempArtifacts(const vector<string>& vTempArtifacts,
                                     const string& sMergedDatabasePath);


int main(int argc, char *argv[])
{
   g_sCometPlusExecutablePath = argv[0];
#ifndef _WIN32
   char szResolvedPath[PATH_MAX];
   if (realpath(argv[0], szResolvedPath) != NULL)
      g_sCometPlusExecutablePath = szResolvedPath;
#endif

   // add git hash to version string if present
   if (strlen(GITHUBSHA) > 0)
   {
      string sTmp = std::string(GITHUBSHA);
      if (sTmp.size() > 7)
         sTmp.resize(7);
      g_sCometVersion = std::string(comet_version) + " (" + sTmp + ")";
   }
   else
      g_sCometVersion = std::string(comet_version);

   if (argc < 2)
      Usage(argv[0]);

   vector<InputFileInfo*> pvInputFiles;
   ICometSearchManager* pCometSearchMgr = GetCometSearchManager();
   char szParamsFile[SIZE_FILE];
   string sMergedDatabasePath;
   vector<string> vTempArtifacts;

   ProcessCmdLine(argc, argv, szParamsFile, pvInputFiles, pCometSearchMgr, sMergedDatabasePath, vTempArtifacts);
   pCometSearchMgr->AddInputFiles(pvInputFiles);

   bool bSearchSucceeded;

   bSearchSucceeded = pCometSearchMgr->DoSearch();

   if (g_bCometPlusKeepTempFiles)
   {
      LogRetainedTempArtifacts(vTempArtifacts, sMergedDatabasePath);
   }
   else
   {
      if (!sMergedDatabasePath.empty())
         remove(sMergedDatabasePath.c_str());

      for (size_t i = 0; i < vTempArtifacts.size(); ++i)
      {
         if (!vTempArtifacts.at(i).empty())
            remove(vTempArtifacts.at(i).c_str());
      }
   }

   CometInterfaces::ReleaseCometSearchManager();

   if (!bSearchSucceeded)
   {
      // We already log errors when search fails, so no need to log the
      // error message again via g_cometStatus
      exit(1);
   }

   return (0);
} // main


void Usage(char *pszCmd,
           bool bFullHelp,
           const char *pszParamsFile)
{
   char szTmp[1024];
   int iSize = sizeof(szTmp);

   logout("\n");
   snprintf(szTmp, iSize, " Comet version \"%s\"\n %s\n", g_sCometVersion.c_str(), copyright);
   logout(szTmp);
   logout("\n");
   snprintf(szTmp, iSize, " CometPlus usage:  %s [options] <input_files>\n", pszCmd);
   logout(szTmp);
   logout("\n");
   logout(" Supported input formats include mzXML, mzML, mzMLb, mzXML.gz, mzML.gz, Thermo raw, mgf, mgf.gz, and ms2 variants (cms2, bms2, ms2)\n");
   logout(" Notes: gzip inputs currently supported are mzXML.gz, mzML.gz, and mgf.gz; ms2.gz is not supported.\n");

   logout("\n");
   logout("       options:  -p         to print out a comet.params.new file\n");
   logout("                 -q         to print out a comet.params.new file with more parameter entries\n");
   logout("                 -P<params> to specify an alternate parameters file (default comet.params)\n");
   logout("                 --params <params>      alias for -P<params>\n");
   logout("                 -N<name>   to specify an alternate output base name; valid for one input, or merged multi-input novel mode\n");
   logout("                 --name <name>          alias for -N<name>\n");
   logout("                 -D<dbase>  to specify a sequence database, overriding entry in parameters file\n");
   logout("                 --database <dbase>     repeatable; supports multi FASTA or multi .idx\n");
   logout("                 --output-folder <dir>  output directory (default: current directory)\n");
   logout("                 --novel_protein <file> novel protein FASTA; digested using Comet settings\n");
   logout("                 --novel_peptide <file> novel peptide input (FASTA or tokenized text)\n");
   logout("                 --output_internal_novel_peptide <file> write internal novel peptide TSV\n");
   logout("                 --internal_novel_peptide <file> reuse internal novel peptide TSV input\n");
   logout("                 --stop-after-saving-novel-peptide stop after writing internal novel peptide TSV\n");
   logout("                 --keep-tmp            keep temporary artifacts on exit for debugging\n");
   logout("                 --scan <file>          scan filter file; delimiters: comma/space/tab/newline\n");
   logout("                 --scan_numbers <list>  explicit scan list, e.g. 1001,1002,1003\n");
   logout("                 note: novel mode requires known DB via --database or params database_name\n");
   logout("                 --thread <num>         override num_threads parameter\n");
   logout("                 -F<num>    to specify the first/start scan to search, overriding entry in parameters file\n");
   logout("                 --first-scan <num>     alias for -F<num>\n");
   logout("                 -L<num>    to specify the last/end scan to search, overriding entry in parameters file\n");
   logout("                 --last-scan <num>      alias for -L<num>\n");
   logout("                            (-F/-L can be set independently; together define a closed range)\n");
   logout("                 -i         create .idx file for fragment ion indexing\n");
   logout("                 -j         create .idx file for peptide indexing\n");
   logout("                 note: --novel_* and --scan* options are not supported with -i/-j\n");
   logout("                 --help-full           print complete CLI parameter-override help\n");
   logout("\n");
   snprintf(szTmp, iSize, "       example:  %s file1.mzXML file2.mzXML\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s -F1000 -L1500 file1.mzXML    <- to search scans 1000 through 1500\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s -PParams.txt *.mzXML         <- use parameters in the file 'Params.txt'\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s --database db1.fasta --database db2.fasta file1.mzML\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s --database db.fasta file1.mzML.gz\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s --database known.fasta --novel_peptide novel.txt --scan_numbers 1001,1002 file1.mzML\n", pszCmd);
   logout(szTmp);
   snprintf(szTmp, iSize, "            or   %s --database known.fasta --novel_protein novel.fasta --scan scan_ids.txt file1.mzML\n", pszCmd);
   logout(szTmp);

   logout("\n");

   if (bFullHelp)
   {
      const char *pszActiveParams = pszParamsFile == NULL ? "comet.params" : pszParamsFile;
      vector<ParamHelpEntry> vHelpEntries;
      string sErrorMsg;

      logout(" Full help - generic parameter overrides\n");
      logout("   This mode prints the detailed override guide and the full parameter-key list.\n");
      logout("   Use --help for short usage only.\n");
      logout("   syntax: --<param_key> <value>\n");
      logout("           --<param_key>=<value>\n");
      logout("   values with spaces must be quoted as one shell argument.\n");
      logout("   values may include inline comments after '#', which are ignored during parsing.\n");
      logout("   variable mod example:\n");
      logout("      --variable_mod01 \"15.994915 M 0 3 -1 0 0  0.0 # Oxidation (M)\"\n");
      logout("   dedicated options must be used for:\n");
      logout("      database_name -> --database\n");
      logout("      num_threads -> --thread\n");
      logout("      scan_range -> --first-scan/--last-scan\n");
      logout("      spectrum_batch_size -> -B\n");
      logout("   note: [COMET_ENZYME_INFO] entries are not command-line overridable keys.\n");

      snprintf(szTmp, iSize, "   overridable keys are loaded from: %s\n", pszActiveParams);
      logout(szTmp);

      if (CollectParamsHelpEntries(pszActiveParams, vHelpEntries, sErrorMsg))
      {
         logout("   key source: selected params file\n");
      }
      else
      {
         vHelpEntries = GetFallbackHelpEntries();
         string sMsg = "   key source: built-in fallback list (params file not readable)\n";
         sMsg += "   reason:" + sErrorMsg;
         logout(sMsg);
      }

      const set<string> &setDedicatedKeys = GetDedicatedOverrideKeys();
      logout("   keys:\n");
      std::sort(vHelpEntries.begin(),
                vHelpEntries.end(),
                [](const ParamHelpEntry &a, const ParamHelpEntry &b)
                {
                   return a.sName < b.sName;
                });
      for (size_t i = 0; i < vHelpEntries.size(); ++i)
      {
         const ParamHelpEntry &entry = vHelpEntries.at(i);
         if (setDedicatedKeys.find(entry.sName) != setDedicatedKeys.end())
            continue;

         snprintf(szTmp, iSize, "      --%s\n", entry.sName.c_str());
         logout(szTmp);

         string sGuide;
         if (!entry.sValueExample.empty())
            sGuide = "example: " + entry.sValueExample;

         if (!entry.sComment.empty())
         {
            if (!sGuide.empty())
               sGuide += "  # " + entry.sComment;
            else
               sGuide = "# " + entry.sComment;
         }

         if (sGuide.empty())
            sGuide = "example: <value from params-file format>";

         snprintf(szTmp, iSize, "         %s\n", sGuide.c_str());
         logout(szTmp);
      }

      logout("\n");
   }

   exit(1);
}


void SetOptions(char *arg,
      char *szParamsFile,
      int *iPrintParams,
      ICometSearchManager *pSearchMgr)
{
   char szTmp[512];
   char szParamStringVal[512];
   int iSize = sizeof(szParamStringVal);

   switch (arg[1])
   {
      case 'D':   // Alternate sequence database.
         strncpy(szTmp, arg+2, 511);
         szTmp[511]='\0';

         if (strlen(szTmp) == 0 )
            logout("Missing text for parameter option -D<database>.  Ignored.\n");
         else
            pSearchMgr->SetParam("database_name", szTmp, szTmp);
         break;
      case 'P':   // Alternate parameters file.
         strncpy(szTmp, arg+2, 511);
         szTmp[511]='\0';

         if (strlen(szTmp) == 0 )
            logout("Missing text for parameter option -P<params>.  Ignored.\n");
         else
            strcpy(szParamsFile, szTmp);
         break;
      case 'N':   // Set basename of output file (for .out, SQT, and pepXML)
         strncpy(szTmp, arg+2, 511);
         szTmp[511]='\0';

         if (strlen(szTmp) == 0 )
            logout("Missing text for parameter option -N<basename>.  Ignored.\n");
         else
            pSearchMgr->SetOutputFileBaseName(szTmp);
         break;
      case 'F':   // first scan
         if (sscanf(arg+2, "%511s", szTmp) == 0 )
            logout("Missing text for parameter option -F<num>.  Ignored.\n");
         else
         {
            IntRange iScanRange;
            pSearchMgr->GetParamValue("scan_range", iScanRange);
            iScanRange.iStart = atoi(szTmp);
            szParamStringVal[0] = '\0';
            snprintf(szParamStringVal, iSize, "%d %d", iScanRange.iStart, iScanRange.iEnd);
            pSearchMgr->SetParam("scan_range", szParamStringVal, iScanRange);
         }
         break;
      case 'L':  // last scan
         if (sscanf(arg+2, "%511s", szTmp) == 0 )
            logout("Missing text for parameter option -L<num>.  Ignored.\n");
         else
         {
            IntRange iScanRange;
            pSearchMgr->GetParamValue("scan_range", iScanRange);
            iScanRange.iEnd = atoi(szTmp);
            szParamStringVal[0] = '\0';
            snprintf(szParamStringVal, iSize, "%d %d", iScanRange.iStart, iScanRange.iEnd);
            pSearchMgr->SetParam("scan_range", szParamStringVal, iScanRange);
         }
         break;
      case 'B':  // batch size
         if (sscanf(arg+2, "%511s", szTmp) == 0 )
            logout("Missing text for parameter option -B<num>.  Ignored.\n");
         else
         {
            int iSpectrumBatchSize = atoi(szTmp);
            szParamStringVal[0] = '\0';
            snprintf(szParamStringVal, iSize, "%d", iSpectrumBatchSize);
            pSearchMgr->SetParam("spectrum_batch_size", szParamStringVal, iSpectrumBatchSize);
         }
         break;
      case 'p':
         *iPrintParams = 1;  // default set of parameters
         break;
      case 'q':
         *iPrintParams = 2;  // include additional parameters such as PEFF, fragment ion index
         break;
      case 'i':
         snprintf(szParamStringVal, iSize, "1");
         pSearchMgr->SetParam("create_fragment_index", szParamStringVal, 1);
         snprintf(szParamStringVal, iSize, "0");
         pSearchMgr->SetParam("create_peptide_index", szParamStringVal, 0);
         break;
      case 'j':
         snprintf(szParamStringVal, iSize, "0");
         pSearchMgr->SetParam("create_fragment_index", szParamStringVal, 0);
         snprintf(szParamStringVal, iSize, "1");
         pSearchMgr->SetParam("create_peptide_index", szParamStringVal, 1);
         break;
      default:
         break;
   }
}


// Parses the command line and determines the type of analysis to perform.
bool ParseCmdLine(char *cmd, InputFileInfo *pInputFile, ICometSearchManager *pSearchMgr)
{
   char *tok;
   char *scan;

   pInputFile->iAnalysisType = 0;

   // Get the file name. Because Windows can have ":" in the file path,
   // we can't just use "strtok" to grab the filename.
   int i;
   int iCmdLen = (int)strlen(cmd);
   for (i = 0; i < iCmdLen; ++i)
   {
      if (cmd[i] == ':')
      {
         if ((i + 1) < iCmdLen)
         {
            if (cmd[i+1] != '\\' && cmd[i+1] != '/')
            {
               break;
            }
         }
      }
   }

   strncpy(pInputFile->szFileName, cmd, i);
   pInputFile->szFileName[i] = '\0';
   if (!ValidateInputFile(pInputFile->szFileName))
   {
      return false;
   }

   // Get additional filters.
   scan = strtok(cmd+i, ":\n");

   // Analyze entire file.
   if (scan == NULL)
   {
      IntRange scanRange;
      if (!pSearchMgr->GetParamValue("scan_range", scanRange))
      {
         scanRange.iStart = 0;
         scanRange.iEnd = 0;
      }

      if (scanRange.iStart == 0 && scanRange.iEnd == 0)
      {
         pInputFile->iAnalysisType = AnalysisType_EntireFile;
         return true;
      }
      else
      {
         pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;

         pInputFile->iFirstScan = scanRange.iStart;
         pInputFile->iLastScan = scanRange.iEnd;

         return true;
      }
   }

   // Analyze a portion of the file.
   if (strchr(scan,'-') != NULL)
   {
      pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;
      tok = strtok(scan, "-\n");
      if (tok != NULL)
         pInputFile->iFirstScan = atoi(tok);
      tok = strtok(NULL,"-\n");
      if (tok != NULL)
         pInputFile->iLastScan = atoi(tok);
   }
   else if (strchr(scan,'+') != NULL)
   {
      pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;
      tok = strtok(scan,"+\n");
      if (tok != NULL)
         pInputFile->iFirstScan = atoi(tok);
      tok = strtok(NULL,"+\n");
      if (tok != NULL)
         pInputFile->iLastScan = pInputFile->iFirstScan + atoi(tok);
   }
   else
   {
      pInputFile->iAnalysisType = AnalysisType_SpecificScan;
      pInputFile->iFirstScan = atoi(scan);
      pInputFile->iLastScan = pInputFile->iFirstScan;
   }

   return true;
} // ParseCmdLine

static bool IsPathSeparatorLocal(char c)
{
   return c == '/' || c == '\\';
}

static bool ContainsPathSeparator(const string& sPath)
{
   for (size_t i = 0; i < sPath.length(); ++i)
   {
      if (IsPathSeparatorLocal(sPath[i]))
         return true;
   }
   return false;
}

static bool IsAbsolutePathLocal(const string& sPath)
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

static string JoinPathLocal(const string& sDir, const string& sName)
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

static string GetDirNameLocal(const string& sPath)
{
   size_t iPos = sPath.find_last_of("/\\");
   if (iPos == string::npos)
      return "";
   if (iPos == 0)
      return sPath.substr(0, 1);
   return sPath.substr(0, iPos);
}

static bool FileExistsLocal(const string& sPath)
{
   struct stat st;
   return stat(sPath.c_str(), &st) == 0;
}

static string NormalizePathKeyLocal(const string& sPath)
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

static bool EndsWithIgnoreCaseLocal(const string& sValue, const string& sSuffix)
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

static bool IsMzMLbPath(const string& sPath)
{
   return EndsWithIgnoreCaseLocal(sPath, ".mzMLb");
}

static bool HasMzMLbInputs(const vector<InputFileInfo*>& vInputs)
{
   for (size_t i = 0; i < vInputs.size(); ++i)
   {
      if (vInputs.at(i) != NULL && IsMzMLbPath(vInputs.at(i)->szFileName))
         return true;
   }

   return false;
}

static string ResolveInternalOutputPath(const string& sPath, const string& sOutputFolder)
{
   if (sPath.empty())
      return sPath;

   if (IsAbsolutePathLocal(sPath) || ContainsPathSeparator(sPath))
      return sPath;

   return JoinPathLocal(sOutputFolder, sPath);
}

static string GetLocalTimestampString()
{
   time_t tNow;
   time(&tNow);
   char szTime[64];
   szTime[0] = '\0';
   strftime(szTime, sizeof(szTime), "%Y-%m-%d %H:%M:%S", localtime(&tNow));
   return szTime;
}

static const char* BoolToOnOffString(bool bValue)
{
   return bValue ? "on" : "off";
}

static void LogRetainedTempArtifacts(const vector<string>& vTempArtifacts,
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

static void LogStageTiming(const string& sStage,
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

static bool PathHasSeparatorForName(const string& sName)
{
   return ContainsPathSeparator(sName);
}

static void AppendPlannedOutputFilesForInput(const InputFileInfo& inputInfo,
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


void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    ICometSearchManager *pSearchMgr,
                    string &sMergedDatabasePath,
                    vector<string> &vTempArtifacts)
{
   const auto tProgramStart = std::chrono::steady_clock::now();

   if (argc == 1)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - no input files specified so nothing to do.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sParamsFile = "comet.params";
   int iPrintParams = 0;
   bool bCreateFragmentIndex = false;
   bool bCreatePeptideIndex = false;
   bool bSetOutputName = false;
   bool bSetFirstScan = false;
   bool bSetLastScan = false;
   bool bSetBatchSize = false;
   bool bSetThreadOverride = false;
   bool bHelpFull = false;
   bool bOutputFolderSpecified = false;
   int iFirstScan = 0;
   int iLastScan = 0;
   int iBatchSize = 0;
   int iThreadOverride = 0;
   string sOutputName;
   string sResolvedOutputInternalNovelPath;
   vector<string> vDatabases;
   vector<string> vInputArgs;
   vector<CmdParamOverride> vCliParamOverrides;
   NovelModeOptions novelOpts;
   g_bCometPlusKeepTempFiles = false;

   for (int iArg = 1; iArg < argc; ++iArg)
   {
      string sArg = argv[iArg];

      if (sArg == "--help" || sArg == "-h")
      {
         Usage(argv[0]);
      }

      if (sArg.rfind("--", 0) == 0)
      {
         string sName;
         string sValue;
         size_t iPos = sArg.find('=');
         if (iPos == string::npos)
         {
            sName = sArg.substr(2);
         }
         else
         {
            sName = sArg.substr(2, iPos - 2);
            sValue = sArg.substr(iPos + 1);
         }

         if (sName.empty())
         {
            string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
            logerr(strErrorMsg);
            exit(1);
         }

         auto RequireValue = [&](const string& sOptName) -> string
         {
            if (!sValue.empty())
               return sValue;

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option --" + sOptName + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         if (sName == "database")
         {
            vDatabases.push_back(RequireValue(sName));
         }
         else if (sName == "help-full")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --help-full does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            bHelpFull = true;
         }
         else if (sName == "thread")
         {
            iThreadOverride = atoi(RequireValue(sName).c_str());
            bSetThreadOverride = true;
         }
         else if (sName == "params")
         {
            sParamsFile = RequireValue(sName);
         }
         else if (sName == "name")
         {
            sOutputName = RequireValue(sName);
            bSetOutputName = true;
         }
         else if (sName == "output-folder")
         {
            novelOpts.sOutputFolder = RequireValue(sName);
            bOutputFolderSpecified = true;
         }
         else if (sName == "first-scan")
         {
            iFirstScan = atoi(RequireValue(sName).c_str());
            bSetFirstScan = true;
         }
         else if (sName == "last-scan")
         {
            iLastScan = atoi(RequireValue(sName).c_str());
            bSetLastScan = true;
         }
         else if (sName == "novel_protein")
         {
            novelOpts.sNovelProteinPath = RequireValue(sName);
         }
         else if (sName == "novel_peptide")
         {
            novelOpts.sNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "output_internal_novel_peptide")
         {
            novelOpts.sOutputInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "internal_novel_peptide")
         {
            novelOpts.sInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "stop-after-saving-novel-peptide")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --stop-after-saving-novel-peptide does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            novelOpts.bStopAfterSavingNovelPeptide = true;
         }
         else if (sName == "keep-tmp")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --keep-tmp does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            g_bCometPlusKeepTempFiles = true;
         }
         else if (sName == "scan")
         {
            novelOpts.sScanFilePath = RequireValue(sName);
         }
         else if (sName == "scan_numbers")
         {
            novelOpts.sScanNumbers = RequireValue(sName);
         }
         else
         {
            CmdParamOverride cliOverride;
            cliOverride.sName = sName;
            cliOverride.sValue = RequireValue(sName);
            vCliParamOverrides.push_back(cliOverride);
         }

         continue;
      }

      if (sArg.length() > 1 && sArg[0] == '-')
      {
         auto RequireShortValue = [&](char cOpt) -> string
         {
            if (sArg.length() > 2)
               return sArg.substr(2);

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option -" + string(1, cOpt) + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         switch (sArg[1])
         {
            case 'p':
               iPrintParams = 1;
               break;
            case 'q':
               iPrintParams = 2;
               break;
            case 'P':
               sParamsFile = RequireShortValue('P');
               break;
            case 'D':
               vDatabases.push_back(RequireShortValue('D'));
               break;
            case 'N':
               sOutputName = RequireShortValue('N');
               bSetOutputName = true;
               break;
            case 'F':
               iFirstScan = atoi(RequireShortValue('F').c_str());
               bSetFirstScan = true;
               break;
            case 'L':
               iLastScan = atoi(RequireShortValue('L').c_str());
               bSetLastScan = true;
               break;
            case 'B':
               iBatchSize = atoi(RequireShortValue('B').c_str());
               bSetBatchSize = true;
               break;
            case 'i':
               bCreateFragmentIndex = true;
               bCreatePeptideIndex = false;
               break;
            case 'j':
               bCreateFragmentIndex = false;
               bCreatePeptideIndex = true;
               break;
            default:
            {
               string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
               logerr(strErrorMsg);
               exit(1);
            }
         }

         continue;
      }

      vInputArgs.push_back(sArg);
   }

   if (bHelpFull)
   {
      Usage(argv[0], true, sParamsFile.c_str());
   }

   if (iPrintParams)
   {
      PrintParams(iPrintParams);
      exit(0);
   }

   strncpy(szParamsFile, sParamsFile.c_str(), SIZE_FILE - 1);
   szParamsFile[SIZE_FILE - 1] = '\0';

   set<string> setParamKeys;
   string sKeyErrMsg;
   if (!CollectParamsFileKeys(szParamsFile, setParamKeys, sKeyErrMsg))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + sKeyErrMsg;
      logerr(strErrorMsg);
      exit(1);
   }

   auto DedicatedHint = [](const string &sParamKey) -> string
   {
      if (sParamKey == "database_name")
         return "--database";
      if (sParamKey == "num_threads")
         return "--thread";
      if (sParamKey == "scan_range")
         return "--first-scan/--last-scan";
      if (sParamKey == "spectrum_batch_size")
         return "-B";
      return "";
   };

   const set<string> &setDedicatedKeys = GetDedicatedOverrideKeys();
   for (size_t i = 0; i < vCliParamOverrides.size(); ++i)
   {
      const string &sOverrideName = vCliParamOverrides.at(i).sName;

      if (setDedicatedKeys.find(sOverrideName) != setDedicatedKeys.end())
      {
         string strErrorMsg = " Error - --" + sOverrideName
            + " is not supported as a generic parameter override; use "
            + DedicatedHint(sOverrideName) + " instead.\n";
         logerr(strErrorMsg);
         exit(1);
      }

         if (setParamKeys.find(sOverrideName) == setParamKeys.end())
         {
            string strErrorMsg = " Error - unknown option \"--" + sOverrideName
               + "\"; key is not present in params file \"" + sParamsFile + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   LoadParameters(szParamsFile, pSearchMgr, vCliParamOverrides);

   if (bSetFirstScan || bSetLastScan)
   {
      IntRange iScanRange;
      pSearchMgr->GetParamValue("scan_range", iScanRange);
      if (bSetFirstScan)
         iScanRange.iStart = iFirstScan;
      if (bSetLastScan)
         iScanRange.iEnd = iLastScan;

      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d %d", iScanRange.iStart, iScanRange.iEnd);
      pSearchMgr->SetParam("scan_range", szParamStringVal, iScanRange);
   }

   if (bSetBatchSize)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iBatchSize);
      pSearchMgr->SetParam("spectrum_batch_size", szParamStringVal, iBatchSize);
   }

   if (bSetThreadOverride)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iThreadOverride);
      pSearchMgr->SetParam("num_threads", szParamStringVal, iThreadOverride);
   }

   if (bCreateFragmentIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "1", 1);
      pSearchMgr->SetParam("create_peptide_index", "0", 0);
   }
   else if (bCreatePeptideIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "0", 0);
      pSearchMgr->SetParam("create_peptide_index", "1", 1);
   }

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty() && !novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --output_internal_novel_peptide requires --novel_protein and/or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.HasInternalNovelInput() && novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --internal_novel_peptide cannot be used together with --novel_protein or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.bStopAfterSavingNovelPeptide && novelOpts.sOutputInternalNovelPeptidePath.empty())
   {
      string strErrorMsg = " Error - --stop-after-saving-novel-peptide requires --output_internal_novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bSetOutputName && bOutputFolderSpecified && PathHasSeparatorForName(sOutputName))
   {
      string strErrorMsg = " Error - when using --name with --output-folder, --name must be a base name without path separators.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if ((bCreateFragmentIndex || bCreatePeptideIndex)
         && (novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()
            || !novelOpts.sOutputInternalNovelPeptidePath.empty()))
   {
      string strErrorMsg = " Error - novel/scan subset options cannot be used with -i or -j.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.sOutputFolder.empty())
      novelOpts.sOutputFolder = ".";

   {
      string sErrorMsg;
      if (!EnsureDirectoryExistsRecursive(novelOpts.sOutputFolder, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   SetCometPlusTempDirectory(novelOpts.sOutputFolder);

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty())
      sResolvedOutputInternalNovelPath = ResolveInternalOutputPath(novelOpts.sOutputInternalNovelPeptidePath,
                                                                   novelOpts.sOutputFolder);

   bool bNeedInputSpectra = !novelOpts.bStopAfterSavingNovelPeptide;
   if ((novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()) && vInputArgs.empty() && bNeedInputSpectra)
   {
      string strErrorMsg = " Error - at least one input spectrum file is required when using novel or scan subset options.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   vector<string> vKnownDatabases = vDatabases;
   if (vKnownDatabases.empty())
   {
      string sDatabaseFromParams;
      if (pSearchMgr->GetParamValue("database_name", sDatabaseFromParams)
            && !sDatabaseFromParams.empty())
      {
         vKnownDatabases.push_back(sDatabaseFromParams);
      }
   }

   bool bKnownAllIdx = true;
   bool bKnownAllFasta = true;
   int iKnownIdxType = 0;
   for (size_t i = 0; i < vKnownDatabases.size(); ++i)
   {
      bool bIdx = IsIdxDatabasePath(vKnownDatabases.at(i));
      if (bIdx)
      {
         bKnownAllFasta = false;
      }
      else
      {
         bKnownAllIdx = false;
      }
   }
   if (!vKnownDatabases.empty() && !(bKnownAllIdx || bKnownAllFasta))
   {
      string strErrorMsg = " Error - mixed database types are not allowed; provide either all FASTA files or all .idx files.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bKnownAllIdx)
   {
      string sProbeError;
      for (size_t i = 0; i < vKnownDatabases.size(); ++i)
      {
         int iObservedType = 0;
         if (!CometPlusProbeIdxType(vKnownDatabases.at(i), iObservedType, sProbeError))
         {
            logerr(sProbeError);
            exit(1);
         }

         if (i == 0)
            iKnownIdxType = iObservedType;
         else if (iObservedType != iKnownIdxType)
         {
            string strErrorMsg = " Error - mixed .idx types are not supported in one search invocation.\n";
            logerr(strErrorMsg);
            exit(1);
         }
      }
   }

   vector<InputFileInfo*> vParsedInputs;
   for (size_t i = 0; i < vInputArgs.size(); ++i)
   {
      InputFileInfo* pInputFileInfo = new InputFileInfo();
      std::vector<char> vInput(vInputArgs.at(i).begin(), vInputArgs.at(i).end());
      vInput.push_back('\0');

      if (!ParseCmdLine(vInput.data(), pInputFileInfo, pSearchMgr))
      {
         string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
            + " Error - input file \"" + std::string(pInputFileInfo->szFileName) + "\" not found.\n";
         logerr(strErrorMsg);
         pvInputFiles.clear();
         exit(1);
      }

      vParsedInputs.push_back(pInputFileInfo);
   }

   const bool bNovelMode = novelOpts.HasNovelMode();
   const bool bUseNovelMergedSearch = (bNovelMode && vParsedInputs.size() > 1);

   if (bSetOutputName && vParsedInputs.size() != 1 && !bUseNovelMergedSearch)
   {
      string strErrorMsg = " Error - --name/-N can be used only when one input spectrum file is provided.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sNovelMergedOutputBase;
   if (bUseNovelMergedSearch)
   {
      string sMergedBaseStem = bSetOutputName ? sOutputName : "cometplus_novel_merged";
      if (sMergedBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive merged output basename for novel multi-input search.\n";
         logerr(strErrorMsg);
         exit(1);
      }

      sNovelMergedOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sMergedBaseStem);
      if (sNovelMergedOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sNovelMergedOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   for (size_t i = 0; i < vParsedInputs.size(); ++i)
   {
      string sInputBaseStem;
      if (bUseNovelMergedSearch)
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);
      else if (bSetOutputName)
         sInputBaseStem = sOutputName;
      else
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);

      if (sInputBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive output basename for input file \"" + string(vParsedInputs.at(i)->szFileName) + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      string sOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sInputBaseStem);
      if (sOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      strncpy(vParsedInputs.at(i)->szBaseName, sOutputBase.c_str(), SIZE_FILE - 1);
      vParsedInputs.at(i)->szBaseName[SIZE_FILE - 1] = '\0';
   }

   int iDecoySearch = 0;
   pSearchMgr->GetParamValue("decoy_search", iDecoySearch);

   int iOutputInt = 0;
   bool bOutputSqtFile = false;
   bool bOutputTxtFile = false;
   bool bOutputPepXmlFile = false;
   int iOutputMzidFile = 0;
   bool bOutputPercolatorFile = false;

   if (pSearchMgr->GetParamValue("output_sqtfile", iOutputInt))
      bOutputSqtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_txtfile", iOutputInt))
      bOutputTxtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_pepxmlfile", iOutputInt))
      bOutputPepXmlFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_mzidentmlfile", iOutputInt))
      iOutputMzidFile = iOutputInt;
   if (pSearchMgr->GetParamValue("output_percolatorfile", iOutputInt))
      bOutputPercolatorFile = (iOutputInt != 0);

   string sOutputSuffix;
   pSearchMgr->GetParamValue("output_suffix", sOutputSuffix);
   string sTxtExt = "txt";
   if (pSearchMgr->GetParamValue("text_file_extension", sTxtExt))
   {
      if (sTxtExt.empty())
         sTxtExt = "txt";
   }

   if (!(bCreateFragmentIndex || bCreatePeptideIndex))
   {
      struct PlannedOutputInput
      {
         InputFileInfo inputInfo;
         string sBaseName;
      };

      vector<PlannedOutputInput> vPlannedOutputInputs;
      if (bUseNovelMergedSearch)
      {
         PlannedOutputInput mergedOutputInput;
         mergedOutputInput.inputInfo.iAnalysisType = AnalysisType_EntireFile;
         mergedOutputInput.sBaseName = sNovelMergedOutputBase;
         vPlannedOutputInputs.push_back(mergedOutputInput);
      }
      else
      {
         for (size_t i = 0; i < vParsedInputs.size(); ++i)
         {
            PlannedOutputInput outputInput;
            outputInput.inputInfo = *vParsedInputs.at(i);
            outputInput.sBaseName = vParsedInputs.at(i)->szBaseName;
            vPlannedOutputInputs.push_back(outputInput);
         }
      }

      vector<string> vPlannedOutputs;
      for (size_t i = 0; i < vPlannedOutputInputs.size(); ++i)
      {
         AppendPlannedOutputFilesForInput(vPlannedOutputInputs.at(i).inputInfo,
                                          vPlannedOutputInputs.at(i).sBaseName,
                                          sOutputSuffix,
                                          sTxtExt,
                                          bOutputSqtFile,
                                          bOutputTxtFile,
                                          bOutputPepXmlFile,
                                          iOutputMzidFile,
                                          bOutputPercolatorFile,
                                          iDecoySearch,
                                          vPlannedOutputs);
      }

      if (!sResolvedOutputInternalNovelPath.empty())
         vPlannedOutputs.push_back(sResolvedOutputInternalNovelPath);

      unordered_map<string, vector<string>> mPathToSources;
      unordered_set<string> setExistingSeen;
      vector<string> vExisting;
      for (size_t i = 0; i < vPlannedOutputs.size(); ++i)
      {
         string sPathKey = NormalizePathKeyLocal(vPlannedOutputs.at(i));
         mPathToSources[sPathKey].push_back(vPlannedOutputs.at(i));
         if (FileExistsLocal(vPlannedOutputs.at(i)) && setExistingSeen.insert(sPathKey).second)
            vExisting.push_back(vPlannedOutputs.at(i));
      }

      vector<string> vInternalConflicts;
      for (auto it = mPathToSources.begin(); it != mPathToSources.end(); ++it)
      {
         if (it->second.size() > 1)
         {
            string sMsg = it->second.at(0) + " (planned " + std::to_string(it->second.size()) + " times)";
            vInternalConflicts.push_back(sMsg);
         }
      }

      if (!vInternalConflicts.empty() || !vExisting.empty())
      {
         string strErrorMsg = " Error - output file conflict(s) detected before search.\n";
         if (!vInternalConflicts.empty())
         {
            strErrorMsg += "  Internal planned conflicts:\n";
            for (size_t i = 0; i < vInternalConflicts.size(); ++i)
               strErrorMsg += "   - " + vInternalConflicts.at(i) + "\n";
         }
         if (!vExisting.empty())
         {
            strErrorMsg += "  Existing files on disk:\n";
            for (size_t i = 0; i < vExisting.size(); ++i)
               strErrorMsg += "   - " + vExisting.at(i) + "\n";
         }
         logerr(strErrorMsg);
         exit(1);
      }
   }

   int iEqualIL = 1;
   pSearchMgr->GetParamValue("equal_I_and_L", iEqualIL);
   bool bTreatSameIL = (iEqualIL != 0);

   vector<double> vNovelMasses;
   string sNovelScoringDatabase;
   vector<NovelPeptideRecord> vNovelRecords;

   if (novelOpts.HasNovelMode())
   {
      if (vKnownDatabases.empty())
      {
         string strErrorMsg = " Error - known database must be provided through --database or database_name in params when using novel options.\n";
         logerr(strErrorMsg);
         exit(1);
      }

      auto AppendUniqueProtein = [](vector<string>& vProteinIds, const string& sProteinId)
      {
         if (sProteinId.empty())
            return;
         for (size_t i = 0; i < vProteinIds.size(); ++i)
         {
            if (vProteinIds.at(i) == sProteinId)
               return;
         }
         vProteinIds.push_back(sProteinId);
      };

      auto BuildNoVarModOverrides = []() -> map<string, string>
      {
         map<string, string> mOverrides;
         mOverrides["max_variable_mods_in_peptide"] = "0";
         mOverrides["require_variable_mod"] = "0";
         for (int iMod = 1; iMod <= 15; ++iMod)
         {
            char szParamName[32];
            snprintf(szParamName, sizeof(szParamName), "variable_mod%02d", iMod);
            mOverrides[szParamName] = "0.0 X 0 3 -1 0 0 0.0";
         }
         return mOverrides;
      };

      if (novelOpts.HasNovelInputOptions())
      {
         string sTmpNoVarModParams;
         {
            string sErrorMsg;
            map<string, string> mNoVarOverrides = BuildNoVarModOverrides();
            if (!BuildTemporaryParamsFile(sParamsFile, mNoVarOverrides, sTmpNoVarModParams, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
            vTempArtifacts.push_back(sTmpNoVarModParams);
         }

         set<string> setKnownPeptides;
         {
            auto tStageStart = std::chrono::steady_clock::now();
            vector<PeptideMassEntry> vKnownEntries;
            string sErrorMsg;

            if (bKnownAllIdx)
            {
               for (size_t iDb = 0; iDb < vKnownDatabases.size(); ++iDb)
               {
                  const string& sDb = vKnownDatabases.at(iDb);

                  int iType = 0;
                  string sProbeError;
                  if (!CometPlusProbeIdxType(sDb, iType, sProbeError))
                  {
                     logerr(sProbeError);
                     exit(1);
                  }

                  vector<PeptideMassEntry> vEntries;
                  bool bOk = false;
                  if (iType == 1)
                     bOk = ParseFragmentIdxEntries(sDb, vEntries, sErrorMsg);
                  else if (iType == 2)
                     bOk = ParsePeptideIdxEntries(sDb, vEntries, false, sErrorMsg);
                  else
                  {
                     sErrorMsg = " Error - unknown .idx type encountered while parsing known database.\n";
                     bOk = false;
                  }

                  if (!bOk)
                  {
                     logerr(sErrorMsg);
                     exit(1);
                  }

                  vKnownEntries.insert(vKnownEntries.end(), vEntries.begin(), vEntries.end());
               }
            }
            else
            {
               string sTmpKnownFasta;
               if (!BuildMergedFasta(vKnownDatabases, sTmpKnownFasta, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }
               vTempArtifacts.push_back(sTmpKnownFasta);

               if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                               sTmpNoVarModParams,
                                               sTmpKnownFasta,
                                               false,
                                               iThreadOverride,
                                               sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }

               string sTmpKnownIdx = sTmpKnownFasta + ".idx";
               vTempArtifacts.push_back(sTmpKnownIdx);
               if (!ParsePeptideIdxEntries(sTmpKnownIdx, vKnownEntries, false, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }
            }

            for (size_t i = 0; i < vKnownEntries.size(); ++i)
            {
               string sKey = NormalizePeptideForCompare(vKnownEntries.at(i).sPeptide, bTreatSameIL);
               if (!sKey.empty())
                  setKnownPeptides.insert(sKey);
            }

            LogStageTiming("known peptide extraction", tStageStart, tProgramStart);
         }

         struct NovelCandidateAggregate
         {
            string sRepresentativePeptide;
            vector<string> vProteinIds;
         };

         unordered_map<string, NovelCandidateAggregate> mNovelCandidates;
         {
            auto tStageStart = std::chrono::steady_clock::now();
            if (!novelOpts.sNovelProteinPath.empty())
            {
               string sErrorMsg;
               string sTmpNovelProtein;
               if (!BuildMergedFasta(vector<string>(1, novelOpts.sNovelProteinPath), sTmpNovelProtein, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }
               vTempArtifacts.push_back(sTmpNovelProtein);

               if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                               sTmpNoVarModParams,
                                               sTmpNovelProtein,
                                               false,
                                               iThreadOverride,
                                               sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }
               string sTmpNovelProteinIdx = sTmpNovelProtein + ".idx";
               vTempArtifacts.push_back(sTmpNovelProteinIdx);

               vector<PeptideMassEntry> vTmp;
               if (!ParsePeptideIdxEntries(sTmpNovelProteinIdx, vTmp, true, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }

               for (size_t i = 0; i < vTmp.size(); ++i)
               {
                  string sNormPep = NormalizePeptideToken(vTmp.at(i).sPeptide);
                  string sKey = NormalizePeptideForCompare(sNormPep, bTreatSameIL);
                  if (sKey.empty())
                     continue;

                  auto it = mNovelCandidates.find(sKey);
                  if (it == mNovelCandidates.end())
                  {
                     NovelCandidateAggregate agg;
                     agg.sRepresentativePeptide = sNormPep;
                     mNovelCandidates[sKey] = agg;
                     it = mNovelCandidates.find(sKey);
                  }

                  for (size_t iProt = 0; iProt < vTmp.at(i).vProteinIds.size(); ++iProt)
                     AppendUniqueProtein(it->second.vProteinIds, vTmp.at(i).vProteinIds.at(iProt));
               }
            }

            if (!novelOpts.sNovelPeptidePath.empty())
            {
               string sErrorMsg;
               vector<string> vInputPeptides;
               if (!ParseNovelPeptideFile(novelOpts.sNovelPeptidePath, vInputPeptides, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }

               if (vInputPeptides.empty())
               {
                  string strErrorMsg = " Error - no peptide entries were parsed from --novel_peptide input.\n";
                  logerr(strErrorMsg);
                  exit(1);
               }

               for (size_t i = 0; i < vInputPeptides.size(); ++i)
               {
                  string sNormPep = NormalizePeptideToken(vInputPeptides.at(i));
                  string sKey = NormalizePeptideForCompare(sNormPep, bTreatSameIL);
                  if (sKey.empty())
                     continue;

                  auto it = mNovelCandidates.find(sKey);
                  if (it == mNovelCandidates.end())
                  {
                     NovelCandidateAggregate agg;
                     agg.sRepresentativePeptide = sNormPep;
                     mNovelCandidates[sKey] = agg;
                     it = mNovelCandidates.find(sKey);
                  }

                  AppendUniqueProtein(it->second.vProteinIds, sNormPep);
               }
            }

            LogStageTiming("novel candidate assembly", tStageStart, tProgramStart);
         }

         {
            auto tStageStart = std::chrono::steady_clock::now();
            int iSubtractedCount = 0;
            vector<string> vKeys;
            vKeys.reserve(mNovelCandidates.size());
            for (auto it = mNovelCandidates.begin(); it != mNovelCandidates.end(); ++it)
               vKeys.push_back(it->first);
            std::sort(vKeys.begin(), vKeys.end());

            vNovelRecords.clear();
            for (size_t i = 0; i < vKeys.size(); ++i)
            {
               if (setKnownPeptides.find(vKeys.at(i)) != setKnownPeptides.end())
               {
                  iSubtractedCount++;
                  continue;
               }

               NovelPeptideRecord rec;
               rec.sPeptide = mNovelCandidates[vKeys.at(i)].sRepresentativePeptide;
               rec.vProteinIds = mNovelCandidates[vKeys.at(i)].vProteinIds;
               vNovelRecords.push_back(rec);
            }

            std::sort(vNovelRecords.begin(),
                      vNovelRecords.end(),
                      [](const NovelPeptideRecord& a, const NovelPeptideRecord& b)
                      {
                         return a.sPeptide < b.sPeptide;
                      });

            for (size_t i = 0; i < vNovelRecords.size(); ++i)
               vNovelRecords.at(i).sPeptideId = "COMETPLUS_NOVEL_" + std::to_string(i + 1);

            char szLogBuf[512];
            snprintf(szLogBuf, sizeof(szLogBuf),
                     " Novel mode: %zu unique novel peptides parsed, %d removed by known-db subtraction, %zu retained.\n",
                     mNovelCandidates.size(),
                     iSubtractedCount,
                     vNovelRecords.size());
            logout(szLogBuf);

            LogStageTiming("novel candidate subtraction", tStageStart, tProgramStart);
         }
      }
      else if (novelOpts.HasInternalNovelInput())
      {
         auto tStageStart = std::chrono::steady_clock::now();
         string sErrorMsg;
         if (!ParseInternalNovelPeptideFile(novelOpts.sInternalNovelPeptidePath, vNovelRecords, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         char szLogBuf[512];
         snprintf(szLogBuf, sizeof(szLogBuf),
                  " Novel mode: imported %zu internal novel peptide records from \"%s\".\n",
                  vNovelRecords.size(),
                  novelOpts.sInternalNovelPeptidePath.c_str());
         logout(szLogBuf);
         LogStageTiming("internal novel peptide import", tStageStart, tProgramStart);
      }

      if (!sResolvedOutputInternalNovelPath.empty())
      {
         auto tStageStart = std::chrono::steady_clock::now();

         string sOutputDir = GetDirNameLocal(sResolvedOutputInternalNovelPath);
         if (!sOutputDir.empty())
         {
            string sErrorMsg;
            if (!EnsureDirectoryExistsRecursive(sOutputDir, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
         }

         std::ofstream outFile(sResolvedOutputInternalNovelPath.c_str(),
                               std::ios::out | std::ios::trunc);
         if (!outFile.good())
         {
            string strErrorMsg = " Error - cannot create internal novel peptide file \"" + sResolvedOutputInternalNovelPath + "\".\n";
            logerr(strErrorMsg);
            exit(1);
         }

         outFile << "peptide\tpeptide_id\tprotein_id\n";
         int iAutoId = 1;
         for (size_t i = 0; i < vNovelRecords.size(); ++i)
         {
            string sPeptide = NormalizePeptideToken(vNovelRecords.at(i).sPeptide);
            if (sPeptide.empty())
               continue;

            string sPeptideId = vNovelRecords.at(i).sPeptideId;
            if (sPeptideId.empty())
            {
               sPeptideId = "COMETPLUS_NOVEL_" + std::to_string(iAutoId);
               iAutoId++;
            }

            string sProteinIdField;
            for (size_t j = 0; j < vNovelRecords.at(i).vProteinIds.size(); ++j)
            {
               if (j > 0)
                  sProteinIdField += ";";
               sProteinIdField += vNovelRecords.at(i).vProteinIds.at(j);
            }

            outFile << sPeptide << "\t" << sPeptideId << "\t" << sProteinIdField << "\n";
         }

         outFile.flush();
         if (!outFile.good())
         {
            string strErrorMsg = " Error - failed writing internal novel peptide file \"" + sResolvedOutputInternalNovelPath + "\".\n";
            logerr(strErrorMsg);
            exit(1);
         }

         char szLogBuf[1024];
         snprintf(szLogBuf,
                  sizeof(szLogBuf),
                  " Saved internal novel peptide file: %s (%zu rows)\n",
                  sResolvedOutputInternalNovelPath.c_str(),
                  vNovelRecords.size());
         logout(szLogBuf);
         LogStageTiming("internal novel peptide export", tStageStart, tProgramStart);

         if (novelOpts.bStopAfterSavingNovelPeptide)
         {
            char szStopBuf[1024];
            snprintf(szStopBuf,
                     sizeof(szStopBuf),
                     " [%s] stop-after-saving-novel-peptide enabled; exiting without spectrum prefilter/search.\n",
                     GetLocalTimestampString().c_str());
            logout(szStopBuf);

            if (g_bCometPlusKeepTempFiles)
            {
               LogRetainedTempArtifacts(vTempArtifacts, sMergedDatabasePath);
            }
            else
            {
               for (size_t i = 0; i < vTempArtifacts.size(); ++i)
               {
                  if (!vTempArtifacts.at(i).empty())
                     remove(vTempArtifacts.at(i).c_str());
               }
               if (!sMergedDatabasePath.empty())
                  remove(sMergedDatabasePath.c_str());
            }
            exit(0);
         }
      }

      if (iDecoySearch == 0)
      {
         logout(" Warning - decoy_search=0 in novel mode; novel peptide inputs are treated as target entries unless already encoded as decoy sequences.\n");
      }

      if (!vNovelRecords.empty())
      {
         auto tStageStart = std::chrono::steady_clock::now();
         string sNovelFastaPath;
         string sErrorMsg;
         if (!WriteNovelRecordsToFasta(vNovelRecords, "cometplus_novel_scoring", sNovelFastaPath, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
         vTempArtifacts.push_back(sNovelFastaPath);

         int iNoCutEnzyme = -1;
         if (!FindNoCutEnzymeNumber(sParamsFile, iNoCutEnzyme, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         map<string, string> mNoCutOverrides;
         mNoCutOverrides["search_enzyme_number"] = std::to_string(iNoCutEnzyme);
         mNoCutOverrides["search_enzyme2_number"] = "0";
         mNoCutOverrides["allowed_missed_cleavage"] = "0";
         mNoCutOverrides["num_enzyme_termini"] = "2";

         string sTmpNoCutParams;
         if (!BuildTemporaryParamsFile(sParamsFile, mNoCutOverrides, sTmpNoCutParams, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
         vTempArtifacts.push_back(sTmpNoCutParams);

         if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                         sTmpNoCutParams,
                                         sNovelFastaPath,
                                         false,
                                         iThreadOverride,
                                         sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         string sTmpNovelMassIdx = sNovelFastaPath + ".idx";
         vTempArtifacts.push_back(sTmpNovelMassIdx);

         vector<PeptideMassEntry> vNovelMassEntries;
         if (!ParsePeptideIdxEntries(sTmpNovelMassIdx, vNovelMassEntries, false, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         set<double> setNovelMasses;
         for (size_t i = 0; i < vNovelMassEntries.size(); ++i)
         {
            if (vNovelMassEntries.at(i).dMass > 0.0)
               setNovelMasses.insert(vNovelMassEntries.at(i).dMass);
         }
         vNovelMasses.assign(setNovelMasses.begin(), setNovelMasses.end());

         if (bKnownAllIdx)
         {
            bool bFragment = (iKnownIdxType == 1);
            if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                            sParamsFile,
                                            sNovelFastaPath,
                                            bFragment,
                                            iThreadOverride,
                                            sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }

            sNovelScoringDatabase = sNovelFastaPath + ".idx";
            vTempArtifacts.push_back(sNovelScoringDatabase);
         }
         else
         {
            sNovelScoringDatabase = sNovelFastaPath;
         }

         LogStageTiming("novel mass calculation", tStageStart, tProgramStart);
      }
      else
      {
         logout(" Warning - no novel peptides remain after subtraction against known database(s); spectrum filtering will likely remove all scans.\n");
      }
   }

   vector<string> vFinalDatabases = vKnownDatabases;
   if (!sNovelScoringDatabase.empty())
      vFinalDatabases.push_back(sNovelScoringDatabase);

   if (!vFinalDatabases.empty())
   {
      string sErrorMsg;
      if (!ConfigureDatabaseInputs(vFinalDatabases, pSearchMgr, sMergedDatabasePath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      if (!sMergedDatabasePath.empty())
         vTempArtifacts.push_back(sMergedDatabasePath);
   }

   char szNovelOutputOnlyParam[32];
   snprintf(szNovelOutputOnlyParam, sizeof(szNovelOutputOnlyParam), "%d", bNovelMode ? 1 : 0);
   pSearchMgr->SetParam("cometplus_novel_output_only", szNovelOutputOnlyParam, bNovelMode ? 1 : 0);
   CometPlusSetNovelOutputOnly(bNovelMode);

   bool bNeedScanPrefilter = bNovelMode || novelOpts.HasExplicitScanFilter();
   if (!bNeedScanPrefilter)
   {
      pvInputFiles = vParsedInputs;

      if (bCreateFragmentIndex || bCreatePeptideIndex)
         return;

      char szSummary[1024];
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu, novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
      logout(szSummary);
      return;
   }

   set<int> setExplicitScans;
   string sErrorMsg;
   if (!novelOpts.sScanNumbers.empty())
   {
      if (!ParseScanIntegersFromString(novelOpts.sScanNumbers, setExplicitScans, sErrorMsg, "--scan_numbers"))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   if (!novelOpts.sScanFilePath.empty())
   {
      if (!ParseScanIntegersFromFile(novelOpts.sScanFilePath, setExplicitScans, sErrorMsg, "--scan"))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   bool bUseExplicitScans = novelOpts.HasExplicitScanFilter();

   NovelMassFilterContext novelMassCtx = {};
   if (!pSearchMgr->GetParamValue("ms_level", novelMassCtx.iMSLevel))
      novelMassCtx.iMSLevel = 2;
   if (bNovelMode)
   {
      if (!BuildNovelMassFilterContext(pSearchMgr, novelMassCtx, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }

   struct PrefilterResult
   {
      string sOriginalPath;
      string sOutputBaseName;
      string sTempMgfPath;
      int iNumScansKept = 0;
      double dElapsedSeconds = 0.0;
      bool bSuccess = false;
      string sErrorMsg;
   };

   int iPrefilterThreads = 1;
   if (bSetThreadOverride && iThreadOverride > 0)
   {
      iPrefilterThreads = iThreadOverride;
   }
   else
   {
      int iNumThreadsParam = 0;
      if (pSearchMgr->GetParamValue("num_threads", iNumThreadsParam) && iNumThreadsParam > 0)
         iPrefilterThreads = iNumThreadsParam;
   }

   if (iPrefilterThreads < 1)
      iPrefilterThreads = 1;
   if (iPrefilterThreads > (int)vParsedInputs.size())
      iPrefilterThreads = (int)vParsedInputs.size();
   if (iPrefilterThreads < 1)
      iPrefilterThreads = 1;

   if (iPrefilterThreads > 1 && HasMzMLbInputs(vParsedInputs))
   {
      iPrefilterThreads = 1;
      char szWarnBuf[1024];
      snprintf(szWarnBuf,
               sizeof(szWarnBuf),
               " [%s] warning: mzMLb input detected; forcing scan prefilter worker_threads=1 for stability.\n",
               GetLocalTimestampString().c_str());
      logout(szWarnBuf);
   }

   char szPrefilterStart[1024];
   snprintf(szPrefilterStart,
            sizeof(szPrefilterStart),
            " [%s] scan prefilter (parallel): %zu input files, worker_threads=%d\n",
            GetLocalTimestampString().c_str(),
            vParsedInputs.size(),
            iPrefilterThreads);
   logout(szPrefilterStart);

   auto tPrefilterStart = std::chrono::steady_clock::now();
   vector<PrefilterResult> vPrefilterResults(vParsedInputs.size());
   std::atomic<size_t> iNextInput(0);
   std::atomic<bool> bPrefilterError(false);
   std::mutex mError;
   string sFirstPrefilterError;

   vector<std::thread> vWorkers;
   vWorkers.reserve((size_t)iPrefilterThreads);
   for (int iWorker = 0; iWorker < iPrefilterThreads; ++iWorker)
   {
      vWorkers.emplace_back([&]()
      {
         while (true)
         {
            size_t iInput = iNextInput.fetch_add(1);
            if (iInput >= vParsedInputs.size())
               break;

            if (bPrefilterError.load())
               break;

            InputFileInfo* pInputFileInfo = vParsedInputs.at(iInput);
            PrefilterResult result;
            result.sOriginalPath = pInputFileInfo->szFileName;
            result.sOutputBaseName = pInputFileInfo->szBaseName;

            auto tInputPrefilterStart = std::chrono::steady_clock::now();
            string sTempMgfPath;
            int iNumScansKept = 0;
            string sLocalError;
            bool bOk = FilterInputFileToTempMgf(*pInputFileInfo,
                                                setExplicitScans,
                                                bUseExplicitScans,
                                                vNovelMasses,
                                                bNovelMode,
                                                novelMassCtx,
                                                sTempMgfPath,
                                                iNumScansKept,
                                                sLocalError);
            result.dElapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - tInputPrefilterStart).count() / 1000.0;

            if (!bOk)
            {
               result.bSuccess = false;
               result.sErrorMsg = sLocalError;
               bPrefilterError.store(true);
               std::lock_guard<std::mutex> lk(mError);
               if (sFirstPrefilterError.empty())
                  sFirstPrefilterError = sLocalError;
            }
            else
            {
               result.sTempMgfPath = sTempMgfPath;
               result.iNumScansKept = iNumScansKept;
               result.bSuccess = true;
            }

            vPrefilterResults.at(iInput) = result;
         }
      });
   }

   for (size_t i = 0; i < vWorkers.size(); ++i)
      vWorkers.at(i).join();

   if (bPrefilterError.load())
   {
      if (sFirstPrefilterError.empty())
      {
         for (size_t i = 0; i < vPrefilterResults.size(); ++i)
         {
            if (!vPrefilterResults.at(i).bSuccess && !vPrefilterResults.at(i).sErrorMsg.empty())
            {
               sFirstPrefilterError = vPrefilterResults.at(i).sErrorMsg;
               break;
            }
         }
      }

      if (sFirstPrefilterError.empty())
         sFirstPrefilterError = " Error - scan prefilter failed.\n";
      logerr(sFirstPrefilterError);
      exit(1);
   }

   for (size_t i = 0; i < vPrefilterResults.size(); ++i)
   {
      vTempArtifacts.push_back(vPrefilterResults.at(i).sTempMgfPath);

      char szLogBuf[1024];
      snprintf(szLogBuf, sizeof(szLogBuf),
               " [%s] scan prefilter: \"%s\" -> %d scans retained (%.3f sec)\n",
               GetLocalTimestampString().c_str(),
               vPrefilterResults.at(i).sOriginalPath.c_str(),
               vPrefilterResults.at(i).iNumScansKept,
               vPrefilterResults.at(i).dElapsedSeconds);
      logout(szLogBuf);
   }

   LogStageTiming("scan prefilter (parallel)", tPrefilterStart, tProgramStart);

   if (bUseNovelMergedSearch)
   {
      auto tMergeStart = std::chrono::steady_clock::now();

      vector<std::pair<string, string>> vMergeInputs;
      vMergeInputs.reserve(vPrefilterResults.size());
      for (size_t i = 0; i < vPrefilterResults.size(); ++i)
      {
         // Keep original TITLE values from filtered MGF shards when basenames are unique.
         vMergeInputs.push_back(std::make_pair(vPrefilterResults.at(i).sTempMgfPath, ""));
      }

      string sMergedMgfPath;
      if (!MergeFilteredMgfFilesWithSourceTag(vMergeInputs, sMergedMgfPath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      vTempArtifacts.push_back(sMergedMgfPath);

      InputFileInfo* pMergedInputFile = new InputFileInfo();
      pMergedInputFile->iAnalysisType = AnalysisType_EntireFile;
      pMergedInputFile->iFirstScan = 0;
      pMergedInputFile->iLastScan = 0;
      strncpy(pMergedInputFile->szFileName, sMergedMgfPath.c_str(), SIZE_FILE - 1);
      pMergedInputFile->szFileName[SIZE_FILE - 1] = '\0';
      strncpy(pMergedInputFile->szBaseName, sNovelMergedOutputBase.c_str(), SIZE_FILE - 1);
      pMergedInputFile->szBaseName[SIZE_FILE - 1] = '\0';

      pvInputFiles.push_back(pMergedInputFile);

      for (size_t i = 0; i < vParsedInputs.size(); ++i)
      {
         delete vParsedInputs.at(i);
         vParsedInputs.at(i) = NULL;
      }
      vParsedInputs.clear();

      char szMergeLog[2048];
      snprintf(szMergeLog,
               sizeof(szMergeLog),
               " [%s] filtered MGF merge: %zu files -> \"%s\"\n",
               GetLocalTimestampString().c_str(),
               vMergeInputs.size(),
               sMergedMgfPath.c_str());
      logout(szMergeLog);

      LogStageTiming("filtered MGF merge", tMergeStart, tProgramStart);
   }
   else
   {
      for (size_t i = 0; i < vParsedInputs.size(); ++i)
      {
         InputFileInfo* pInputFileInfo = vParsedInputs.at(i);
         strncpy(pInputFileInfo->szFileName, vPrefilterResults.at(i).sTempMgfPath.c_str(), SIZE_FILE - 1);
         pInputFileInfo->szFileName[SIZE_FILE - 1] = '\0';
         strncpy(pInputFileInfo->szBaseName, vPrefilterResults.at(i).sOutputBaseName.c_str(), SIZE_FILE - 1);
         pInputFileInfo->szBaseName[SIZE_FILE - 1] = '\0';
         pvInputFiles.push_back(pInputFileInfo);
      }
   }

   char szSummary[1024];
   if (bUseNovelMergedSearch)
   {
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu (merged from %zu source files), novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               vPrefilterResults.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
   }
   else
   {
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu, novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
   }
   logout(szSummary);
} // ProcessCmdLine


void PrintParams(int iPrintParams)
{
   FILE *fp;

   if ( (fp=fopen("comet.params.new", "w"))==NULL)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - cannot write file comet.params.new\n";
      logerr(strErrorMsg);
      exit(1);
   }

   fprintf(fp, "# comet_version %s\n\
# Comet MS/MS search engine parameters file.\n\
# Everything following the '#' symbol is treated as a comment.\n", g_sCometVersion.c_str());

   fprintf(fp,
"#\n\
database_name = /some/path/db.fasta\n\
decoy_search = 0                       # 0=no (default), 1=internal decoy concatenated, 2=internal decoy separate\n\
\n\
num_threads = 0                        # 0=poll CPU to set num threads; else specify num threads directly (max %d)\n\n", MAX_THREADS);

   if (iPrintParams == 2)
   {
      fprintf(fp, "\nspectral_library_name = /some/path/speclib.file\n");
      fprintf(fp, "spectral_library_ms_level = 1\n\n");

      fprintf(fp,
"#\n\
# PEFF - PSI Extended FASTA Format\n\
#\n\
peff_format = 0                        # 0=no (normal fasta, default), 1=PEFF PSI-MOD, 2=PEFF Unimod\n\
peff_obo =                             # path to PSI Mod or Unimod OBO file\n\
\n\
#\n\
# fragment ion index; limited to 5 variable mods and up to 5 modified residues per mod\n\
#\n\
fragindex_min_ions_score = 3           # minimum number of matched fragment ion index peaks for scoring\n\
fragindex_min_ions_report = 3          # minimum number of matched fragment ion index peaks for reporting(>= fragindex_min_ions_score)\n\
fragindex_num_spectrumpeaks = 100      # number of peaks from spectrum to use for fragment ion index matching\n\
fragindex_min_fragmentmass = 200.0     # low mass cutoff for fragment ions\n\
fragindex_max_fragmentmass = 2000.0    # high mass cutoff for fragment ions\n\
fragindex_skipreadprecursors = 0       # 0=read precursors to limit fragment ion index, 1=skip reading precursors\n\n");
   }

   fprintf(fp,
"#\n\
# masses\n\
#\n\
peptide_mass_tolerance_upper = 20.0    # upper bound of the precursor mass tolerance\n\
peptide_mass_tolerance_lower = -20.0   # lower bound of the precursor mass tolerance; USUALLY NEGATIVE TO BE LOWER THAN 0\n\
peptide_mass_units = 2                 # 0=amu, 1=mmu, 2=ppm\n\
precursor_tolerance_type = 1           # 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances\n\
isotope_error = 2                      # 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1\n");

   if (iPrintParams == 2)
   {
      fprintf(fp,
"mass_type_parent = 1                   # 0=average masses, 1=monoisotopic masses\n\
mass_type_fragment = 1                 # 0=average masses, 1=monoisotopic masses\n");
   }

   fprintf(fp,
"\n\
#\n\
# search enzyme\n\
#\n\
search_enzyme_number = 1               # choose from list at end of this params file\n\
search_enzyme2_number = 0              # second enzyme; set to 0 if no second enzyme\n\
sample_enzyme_number = 1               # specifies the sample enzyme which is possibly different than the one applied to the search;\n\
                                       # used by PeptideProphet to calculate NTT & NMC in pepXML output (default=1 for trypsin).\n\
num_enzyme_termini = 2                 # 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific\n\
allowed_missed_cleavage = 2            # maximum value is 5; for enzyme search\n\
\n\
#\n\
# Up to 15 variable_mod entries are supported for a standard search; manually add additional entries as needed\n\
# format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required> <neutral_loss>\n\
#     e.g. 79.966331 STY 0 3 -1 0 0 97.976896\n\
#\n\
variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0\n\
variable_mod02 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod03 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod04 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod05 = 0.0 X 0 3 -1 0 0 0.0\n");
   if (iPrintParams == 2)
   {
      fprintf(fp,
"variable_mod06 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod07 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod08 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod09 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod10 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod11 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod12 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod13 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod14 = 0.0 X 0 3 -1 0 0 0.0\n\
variable_mod15 = 0.0 X 0 3 -1 0 0 0.0\n");
   }

   fprintf(fp,
"max_variable_mods_in_peptide = 5\n\
require_variable_mod = 0\n");
   if (iPrintParams == 2)
   {
      fprintf(fp, "protein_modslist_file =                # limit variable mods to subset of specified proteins if this file is specified & present\n");
   }

   fprintf(fp,
"\n\
#\n\
# fragment ions\n\
#\n\
# ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1\n\
# high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0, spectrum_batch_size = 15000\n\
#\n\
fragment_bin_tol = 0.02                # binning to use on fragment ions\n\
fragment_bin_offset = 0.0              # offset position to start the binning (0.0 to 1.0)\n\
theoretical_fragment_ions = 0          # 0=use flanking peaks, 1=M peak only\n\
use_A_ions = 0\n\
use_B_ions = 1\n\
use_C_ions = 0\n\
use_X_ions = 0\n\
use_Y_ions = 1\n\
use_Z_ions = 0\n\
use_Z1_ions = 0\n\
use_NL_ions = 0                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks\n\
\n\
#\n\
# output\n\
#\n\
output_sqtfile = 0                     # 0=no, 1=yes  write sqt file\n\
output_txtfile = 0                     # 0=no, 1=yes, 2=Crux-formatted  write tab-delimited txt file\n\
output_pepxmlfile = 1                  # 0=no, 1=yes  write pepXML file\n\
output_mzidentmlfile = 0               # 0=no, 1=yes  write mzIdentML file\n\
output_percolatorfile = 0              # 0=no, 1=yes  write Percolator pin file\n");

   if (iPrintParams == 2)
   {
      fprintf(fp,
"print_expect_score = 1                 # 0=no, 1=yes to replace Sp with expect in out & sqt\n\
print_ascorepro_score = 1              # 0=no, 0 to 5 to localize variable_mod01 to _mod05; -1 to localize all variable mods\n");
   }
 
   fprintf(fp,
"num_output_lines = 5                   # num peptide results to show\n\
\n\
#\n\
# mzXML/mzML/raw file parameters\n\
#\n\
scan_range = 0 0                       # start and end scan range to search; either entry can be set independently\n\
precursor_charge = 0 0                 # precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter\n\
override_charge = 0                    # 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online\n\
ms_level = 2                           # MS level to analyze, valid are levels 2 (default) or 3\n\
activation_method = ALL                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, ETD+SA, PQD, HCD, IRMPD, SID\n\
\n\
#\n\
# misc parameters\n\
#\n\
digest_mass_range = 600.0 5000.0       # MH+ peptide mass range to analyze\n\
peptide_length_range = 5 50            # minimum and maximum peptide length to analyze (default min 1 to allowed max %d)\n",
      MAX_PEPTIDE_LEN);

   if (iPrintParams == 2)
   {
      fprintf(fp,
"pinfile_protein_delimiter =            # blank = default 'tab' delimiter between proteins; enter a char/string to use in place of the tab; Percolator pin output only\n");
      fprintf(fp,
"num_results = 100                      # number of results to store internally for Sp rank only; if Sp rank is not used, set this to num_output_lines\n");
   }

fprintf(fp,
"max_duplicate_proteins = 10            # maximum number of additional duplicate protein names to report for each peptide ID; -1 reports all duplicates\n\
max_fragment_charge = 3                # set maximum fragment charge state to analyze (allowed max %d)\n\
min_precursor_charge = 1               # set minimum precursor charge state to analyze (1 if missing)\n\
max_precursor_charge = 6               # set maximum precursor charge state to analyze (allowed max %d)\n",
      MAX_FRAGMENT_CHARGE,
      MAX_PRECURSOR_CHARGE);

fprintf(fp,
"clip_nterm_methionine = 0              # 0=leave protein sequences as-is; 1=also consider sequence w/o N-term methionine\n\
spectrum_batch_size = 15000            # max. # of spectra to search at a time; 0 to search the entire scan range in one loop\n\
decoy_prefix = DECOY_                  # decoy entries are denoted by this string which is pre-pended to each protein accession\n\
equal_I_and_L = 1                      # 0=treat I and L as different; 1=treat I and L as same\n\
mass_offsets =                         # one or more mass offsets to search (values substracted from deconvoluted precursor mass)\n\
\n\
#\n\
# spectral processing\n\
#\n\
minimum_peaks = 10                     # required minimum number of peaks in spectrum to search (default 10)\n");

fprintf(fp,
"minimum_intensity = 0                 # minimum intensity value to read in\n\
remove_precursor_peak = 0              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD), 3=phosphate neutral loss peaks\n\
remove_precursor_tolerance = 1.5       # +- Da tolerance for precursor removal\n\
clear_mz_range = 0.0 0.0               # clear out all peaks in the specified m/z range e.g. remove reporter ion region of TMT spectra\n\
percentage_base_peak = 0.0             # specify a percentage (e.g. \"0.05\" for 5%%) of the base peak intensity as a minimum intensity threshold\n\
\n\
#\n\
# static modifications\n\
#\n\
add_Cterm_peptide = 0.0\n\
add_Nterm_peptide = 0.0\n\
add_Cterm_protein = 0.0\n\
add_Nterm_protein = 0.0\n\
\n\
add_G_glycine = 0.0000                 # added to G - avg.  57.0513, mono.  57.02146\n\
add_A_alanine = 0.0000                 # added to A - avg.  71.0779, mono.  71.03711\n\
add_S_serine = 0.0000                  # added to S - avg.  87.0773, mono.  87.03203\n\
add_P_proline = 0.0000                 # added to P - avg.  97.1152, mono.  97.05276\n\
add_V_valine = 0.0000                  # added to V - avg.  99.1311, mono.  99.06841\n\
add_T_threonine = 0.0000               # added to T - avg. 101.1038, mono. 101.04768\n\
add_C_cysteine = 57.021464             # added to C - avg. 103.1429, mono. 103.00918\n\
add_L_leucine = 0.0000                 # added to L - avg. 113.1576, mono. 113.08406\n\
add_I_isoleucine = 0.0000              # added to I - avg. 113.1576, mono. 113.08406\n\
add_N_asparagine = 0.0000              # added to N - avg. 114.1026, mono. 114.04293\n\
add_D_aspartic_acid = 0.0000           # added to D - avg. 115.0874, mono. 115.02694\n\
add_Q_glutamine = 0.0000               # added to Q - avg. 128.1292, mono. 128.05858\n\
add_K_lysine = 0.0000                  # added to K - avg. 128.1723, mono. 128.09496\n\
add_E_glutamic_acid = 0.0000           # added to E - avg. 129.1140, mono. 129.04259\n\
add_M_methionine = 0.0000              # added to M - avg. 131.1961, mono. 131.04048\n\
add_H_histidine = 0.0000               # added to H - avg. 137.1393, mono. 137.05891\n\
add_F_phenylalanine = 0.0000           # added to F - avg. 147.1739, mono. 147.06841\n\
add_U_selenocysteine = 0.0000          # added to U - avg. 150.0379, mono. 150.95363\n\
add_R_arginine = 0.0000                # added to R - avg. 156.1857, mono. 156.10111\n\
add_Y_tyrosine = 0.0000                # added to Y - avg. 163.0633, mono. 163.06333\n\
add_W_tryptophan = 0.0000              # added to W - avg. 186.0793, mono. 186.07931\n\
add_O_pyrrolysine = 0.0000             # added to O - avg. 237.2982, mono  237.14773\n\
add_B_user_amino_acid = 0.0000         # added to B - avg.   0.0000, mono.   0.00000\n\
add_J_user_amino_acid = 0.0000         # added to J - avg.   0.0000, mono.   0.00000\n\
add_X_user_amino_acid = 0.0000         # added to X - avg.   0.0000, mono.   0.00000\n\
add_Z_user_amino_acid = 0.0000         # added to Z - avg.   0.0000, mono.   0.00000\n\n");

   if (0)  // do not print these parameters out
   {
fprintf(fp,
"#\n\
# These set_X_residue parameters will override the default AA masses for both precursor and fragment calculations.\n\
# They are applied if the parameter value is not zero.\n\
#\n\
set_G_glycine = 0.0000\n\
set_A_alanine = 0.0000\n\
set_S_serine = 0.0000\n\
set_P_proline = 0.0000\n\
set_V_valine = 0.0000\n\
set_T_threonine = 0.0000\n\
set_C_cysteine = 0.0000\n\
set_L_leucine = 0.0000\n\
set_I_isoleucine = 0.0000\n\
set_N_asparagine = 0.0000\n\
set_D_aspartic_acid = 0.0000\n\
set_Q_glutamine = 0.0000\n\
set_K_lysine = 0.0000\n\
set_E_glutamic_acid = 0.0000\n\
set_M_methionine = 0.0000\n\
set_H_histidine = 0.0000\n\
set_F_phenylalanine = 0.0000\n\
set_U_selenocysteine = 0.0000\n\
set_R_arginine = 0.0000\n\
set_Y_tyrosine = 0.0000\n\
set_W_tryptophan = 0.0000\n\
set_O_pyrrolysine = 0.0000\n\
set_B_user_amino_acid = 0.0000\n\
set_J_user_amino_acid = 0.0000\n\
set_X_user_amino_acid = 0.0000\n\
set_Z_user_amino_acid = 0.0000\n\n");
   }

fprintf(fp,
"#\n\
# COMET_ENZYME_INFO _must_ be at the end of this parameters file\n\
# Enzyme entries can be added/deleted/edited\n\
#\n\
[COMET_ENZYME_INFO]\n\
0.  Cut_everywhere         0      -           -\n\
1.  Trypsin                1      KR          P\n\
2.  Trypsin/P              1      KR          -\n\
3.  Lys_C                  1      K           P\n\
4.  Lys_N                  0      K           -\n\
5.  Arg_C                  1      R           P\n\
6.  Asp_N                  0      DN          -\n\
7.  CNBr                   1      M           -\n\
8.  Asp-N_ambic            1      DE          -\n\
9.  PepsinA                1      FL          -\n\
10. Chymotrypsin           1      FWYL        P\n\
11. No_cut                 1      @           @\n\
\n");

   std::string sTmp = " Comet version \"" + g_sCometVersion + "\"\n " + copyright + "\n\n Created:  comet.params.new\n\n";
   logout(sTmp.c_str());
   fclose(fp);

} // PrintParams


bool ValidateInputFile(char *pszInputFileName)
{
   FILE *fp;
   if ((fp = fopen(pszInputFileName, "r")) == NULL)
   {
      return false;
   }
   fclose(fp);
   return true;
}
