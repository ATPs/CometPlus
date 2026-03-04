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


#include "CometPlusApp.h"
#include "CometPlusRuntimeUtils.h"
#include "CometData.h"
#include "CometDataInternal.h"
#include "CometPlusParams.h"
#include "githubsha.h"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifndef _WIN32
#include <unistd.h>
#endif

using namespace CometInterfaces;

string g_sCometPlusExecutablePath;
bool g_bCometPlusKeepTempFiles = false;

static string ResolveExecutablePathLocal(const char* pszArgv0)
{
   if (pszArgv0 == NULL || *pszArgv0 == '\0')
      return "";

#ifdef _WIN32
   return string(pszArgv0);
#else
#ifdef __linux__
   char szProcExe[PATH_MAX + 1];
   ssize_t iRead = readlink("/proc/self/exe", szProcExe, PATH_MAX);
   if (iRead > 0)
   {
      szProcExe[iRead] = '\0';
      return string(szProcExe);
   }
#endif

   char szResolvedPath[PATH_MAX];
   if (realpath(pszArgv0, szResolvedPath) != NULL)
      return string(szResolvedPath);

   return string(pszArgv0);
#endif
}


int main(int argc, char *argv[])
{
   g_sCometPlusExecutablePath = ResolveExecutablePathLocal(argv[0]);

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
   logout("                 env COMETPLUS_PREFILTER_WORKER can override worker path for mzMLb process-prefilter\n");
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
