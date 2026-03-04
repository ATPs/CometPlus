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
#include "NovelModeUtils.h"
#include <climits>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

using namespace std;

static void PrintUsage(const char* pszCmd)
{
   const char* pszProg = (pszCmd != NULL && *pszCmd != '\0')
      ? pszCmd
      : "cometplus_prefilter_worker";

   fprintf(stdout, "\n");
   fprintf(stdout, "CometPlus prefilter worker\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Purpose:\n");
   fprintf(stdout, "  Execute one spectrum prefilter job and write a machine-readable result file.\n");
   fprintf(stdout, "  This binary is normally launched by cometplus for mzMLb process-parallel prefilter.\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Usage:\n");
   fprintf(stdout, "  %s --job <job.tsv>\n", pszProg);
   fprintf(stdout, "  %s --job=<job.tsv>\n", pszProg);
   fprintf(stdout, "  %s --help\n", pszProg);
   fprintf(stdout, "\n");
   fprintf(stdout, "Options:\n");
   fprintf(stdout, "  --job <path>    required. Worker job file path.\n");
   fprintf(stdout, "  --help, -h      print this help text.\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Job File Format:\n");
   fprintf(stdout, "  TAB-separated key/value text file, one key-value per line:\n");
   fprintf(stdout, "    <key>\\t<value>\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Required keys:\n");
   fprintf(stdout, "  result_file, input_file, temp_dir\n");
   fprintf(stdout, "  analysis_type, first_scan, last_scan\n");
   fprintf(stdout, "  use_explicit_scans, explicit_scans_file\n");
   fprintf(stdout, "  use_novel_mass_filter, novel_masses_file, mass_offsets_file\n");
   fprintf(stdout, "  ctx_dPepMassLow, ctx_dPepMassHigh, ctx_dTolLower, ctx_dTolUpper\n");
   fprintf(stdout, "  ctx_iTolUnits, ctx_iTolType, ctx_iIsotopeError\n");
   fprintf(stdout, "  ctx_iStartCharge, ctx_iEndCharge, ctx_iOverrideCharge\n");
   fprintf(stdout, "  ctx_iMinPrecursorCharge, ctx_iMaxPrecursorCharge\n");
   fprintf(stdout, "  ctx_bCorrectMass, ctx_iMSLevel\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Optional keys:\n");
   fprintf(stdout, "  internal_novel_peptide_file (preferred novel-mass source when detailed TSV columns exist)\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Boolean fields use 0/1.\n");
   fprintf(stdout, "When a filter is disabled, its companion list file key can be empty.\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Result File Format (written to result_file):\n");
   fprintf(stdout, "  success\\t0|1\n");
   fprintf(stdout, "  temp_mgf_path\\t<path>\n");
   fprintf(stdout, "  num_scans_kept\\t<int>\n");
   fprintf(stdout, "  error\\t<single-line error message>\n");
   fprintf(stdout, "\n");
   fprintf(stdout, "Example:\n");
   fprintf(stdout, "  %s --job /tmp/cometplus_prefilter_job.tsv\n", pszProg);
   fprintf(stdout, "\n");
}

static void SanitizeInline(string& sValue)
{
   for (size_t i = 0; i < sValue.length(); ++i)
   {
      if (sValue[i] == '\n' || sValue[i] == '\r' || sValue[i] == '\t')
         sValue[i] = ' ';
   }
}

static bool ParseIntStrict(const string& sValue, int& iOut)
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

static bool ParseDoubleStrict(const string& sValue, double& dOut)
{
   char* pEnd = NULL;
   errno = 0;
   double dValue = strtod(sValue.c_str(), &pEnd);
   if (pEnd == NULL || *pEnd != '\0' || errno == ERANGE)
      return false;
   dOut = dValue;
   return true;
}

static bool ParseBool01(const string& sValue, bool& bOut)
{
   int iValue = 0;
   if (!ParseIntStrict(sValue, iValue))
      return false;
   if (iValue != 0 && iValue != 1)
      return false;
   bOut = (iValue != 0);
   return true;
}

static bool ReadDoubleListFile(const string& sPath,
                               vector<double>& vValues,
                               string& sErrorMsg)
{
   vValues.clear();
   sErrorMsg.clear();

   if (sPath.empty())
      return true;

   ifstream inFile(sPath.c_str(), ios::in | ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read worker list file \"" + sPath + "\".\n";
      return false;
   }

   string sLine;
   while (std::getline(inFile, sLine))
   {
      if (!sLine.empty() && sLine[sLine.length() - 1] == '\r')
         sLine.erase(sLine.length() - 1);
      if (sLine.empty())
         continue;

      double dValue = 0.0;
      if (!ParseDoubleStrict(sLine, dValue))
      {
         sErrorMsg = " Error - invalid floating-point token \"" + sLine + "\" in \"" + sPath + "\".\n";
         return false;
      }
      vValues.push_back(dValue);
   }

   return true;
}

static bool ReadJobFile(const string& sPath,
                        map<string, string>& mValues,
                        string& sErrorMsg)
{
   mValues.clear();
   sErrorMsg.clear();

   ifstream inFile(sPath.c_str(), ios::in | ios::binary);
   if (!inFile.good())
   {
      sErrorMsg = " Error - cannot read prefilter worker job file \"" + sPath + "\".\n";
      return false;
   }

   string sLine;
   while (std::getline(inFile, sLine))
   {
      if (!sLine.empty() && sLine[sLine.length() - 1] == '\r')
         sLine.erase(sLine.length() - 1);
      if (sLine.empty())
         continue;

      size_t iTab = sLine.find('\t');
      if (iTab == string::npos)
      {
         sErrorMsg = " Error - malformed worker job line: \"" + sLine + "\".\n";
         return false;
      }

      string sKey = sLine.substr(0, iTab);
      string sValue = sLine.substr(iTab + 1);
      mValues[sKey] = sValue;
   }

   return true;
}

static bool GetRequiredField(const map<string, string>& mValues,
                             const string& sKey,
                             string& sOutValue,
                             string& sErrorMsg)
{
   auto it = mValues.find(sKey);
   if (it == mValues.end())
   {
      sErrorMsg = " Error - missing field \"" + sKey + "\" in worker job file.\n";
      return false;
   }

   sOutValue = it->second;
   return true;
}

static bool WriteWorkerResultFile(const string& sResultPath,
                                  bool bSuccess,
                                  const string& sTempMgfPath,
                                  int iNumScansKept,
                                  const string& sErrorMsg)
{
   ofstream outFile(sResultPath.c_str(), ios::out | ios::trunc | ios::binary);
   if (!outFile.good())
      return false;

   outFile << "success\t" << (bSuccess ? "1" : "0") << "\n";
   outFile << "temp_mgf_path\t" << sTempMgfPath << "\n";
   outFile << "num_scans_kept\t" << iNumScansKept << "\n";

   string sCleanError = sErrorMsg;
   SanitizeInline(sCleanError);
   outFile << "error\t" << sCleanError << "\n";

   outFile.flush();
   return outFile.good();
}

static int ExitWithResult(const string& sResultPath,
                          const string& sErrorMsg,
                          int iCode)
{
   string sCleanError = sErrorMsg;
   if (sCleanError.empty())
      sCleanError = " Error - prefilter worker failed.\n";

   if (!sResultPath.empty())
   {
      if (!WriteWorkerResultFile(sResultPath, false, "", 0, sCleanError))
      {
         string sInline = sCleanError;
         SanitizeInline(sInline);
         fprintf(stderr, "%s\n", sInline.c_str());
      }
   }
   else
   {
      string sInline = sCleanError;
      SanitizeInline(sInline);
      fprintf(stderr, "%s\n", sInline.c_str());
   }

   return iCode;
}

int main(int argc, char* argv[])
{
   string sJobPath;

   for (int i = 1; i < argc; ++i)
   {
      string sArg = argv[i];
      if (sArg == "--help" || sArg == "-h")
      {
         PrintUsage(argv[0]);
         return 0;
      }
      else if (sArg == "--job")
      {
         if (i + 1 >= argc)
         {
            fprintf(stderr, "Error - --job requires a value.\n");
            return 1;
         }
         sJobPath = argv[++i];
      }
      else if (sArg.rfind("--job=", 0) == 0)
      {
         sJobPath = sArg.substr(6);
      }
      else
      {
         fprintf(stderr, "Error - unknown option \"%s\".\n", sArg.c_str());
         fprintf(stderr, "Run \"%s --help\" for usage details.\n", argv[0]);
         return 1;
      }
   }

   if (sJobPath.empty())
   {
      fprintf(stderr, "Error - missing required option --job <path>.\n");
      fprintf(stderr, "Run \"%s --help\" for usage details.\n", argv[0]);
      return 1;
   }

   map<string, string> mJob;
   string sErrorMsg;
   if (!ReadJobFile(sJobPath, mJob, sErrorMsg))
   {
      fprintf(stderr, "%s", sErrorMsg.c_str());
      return 1;
   }

   string sResultPath;
   if (!GetRequiredField(mJob, "result_file", sResultPath, sErrorMsg))
   {
      fprintf(stderr, "%s", sErrorMsg.c_str());
      return 1;
   }

   string sInputPath;
   if (!GetRequiredField(mJob, "input_file", sInputPath, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);

   string sTempDir;
   if (!GetRequiredField(mJob, "temp_dir", sTempDir, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);
   SetCometPlusTempDirectory(sTempDir);

   string sTmpValue;
   int iAnalysisType = 0;
   int iFirstScan = 0;
   int iLastScan = 0;
   bool bUseExplicitScans = false;
   bool bUseNovelMassFilter = false;

   if (!GetRequiredField(mJob, "analysis_type", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, iAnalysisType))
      return ExitWithResult(sResultPath, " Error - invalid analysis_type in worker job.\n", 1);

   if (!GetRequiredField(mJob, "first_scan", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, iFirstScan))
      return ExitWithResult(sResultPath, " Error - invalid first_scan in worker job.\n", 1);

   if (!GetRequiredField(mJob, "last_scan", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, iLastScan))
      return ExitWithResult(sResultPath, " Error - invalid last_scan in worker job.\n", 1);

   if (!GetRequiredField(mJob, "use_explicit_scans", sTmpValue, sErrorMsg)
         || !ParseBool01(sTmpValue, bUseExplicitScans))
      return ExitWithResult(sResultPath, " Error - invalid use_explicit_scans in worker job.\n", 1);

   if (!GetRequiredField(mJob, "use_novel_mass_filter", sTmpValue, sErrorMsg)
         || !ParseBool01(sTmpValue, bUseNovelMassFilter))
      return ExitWithResult(sResultPath, " Error - invalid use_novel_mass_filter in worker job.\n", 1);

   string sExplicitScansFile;
   if (!GetRequiredField(mJob, "explicit_scans_file", sExplicitScansFile, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);

   string sNovelMassesFile;
   if (!GetRequiredField(mJob, "novel_masses_file", sNovelMassesFile, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);

   string sInternalNovelPeptideFile;
   auto itInternalNovel = mJob.find("internal_novel_peptide_file");
   if (itInternalNovel != mJob.end())
      sInternalNovelPeptideFile = itInternalNovel->second;

   string sMassOffsetsFile;
   if (!GetRequiredField(mJob, "mass_offsets_file", sMassOffsetsFile, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);

   NovelMassFilterContext ctx = {};

   if (!GetRequiredField(mJob, "ctx_dPepMassLow", sTmpValue, sErrorMsg)
         || !ParseDoubleStrict(sTmpValue, ctx.dPepMassLow))
      return ExitWithResult(sResultPath, " Error - invalid ctx_dPepMassLow in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_dPepMassHigh", sTmpValue, sErrorMsg)
         || !ParseDoubleStrict(sTmpValue, ctx.dPepMassHigh))
      return ExitWithResult(sResultPath, " Error - invalid ctx_dPepMassHigh in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_dTolLower", sTmpValue, sErrorMsg)
         || !ParseDoubleStrict(sTmpValue, ctx.dTolLower))
      return ExitWithResult(sResultPath, " Error - invalid ctx_dTolLower in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_dTolUpper", sTmpValue, sErrorMsg)
         || !ParseDoubleStrict(sTmpValue, ctx.dTolUpper))
      return ExitWithResult(sResultPath, " Error - invalid ctx_dTolUpper in worker job.\n", 1);

   if (!GetRequiredField(mJob, "ctx_iTolUnits", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iTolUnits))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iTolUnits in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iTolType", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iTolType))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iTolType in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iIsotopeError", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iIsotopeError))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iIsotopeError in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iStartCharge", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iStartCharge))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iStartCharge in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iEndCharge", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iEndCharge))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iEndCharge in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iOverrideCharge", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iOverrideCharge))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iOverrideCharge in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iMinPrecursorCharge", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iMinPrecursorCharge))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iMinPrecursorCharge in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iMaxPrecursorCharge", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iMaxPrecursorCharge))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iMaxPrecursorCharge in worker job.\n", 1);

   if (!GetRequiredField(mJob, "ctx_bCorrectMass", sTmpValue, sErrorMsg)
         || !ParseBool01(sTmpValue, ctx.bCorrectMass))
      return ExitWithResult(sResultPath, " Error - invalid ctx_bCorrectMass in worker job.\n", 1);
   if (!GetRequiredField(mJob, "ctx_iMSLevel", sTmpValue, sErrorMsg)
         || !ParseIntStrict(sTmpValue, ctx.iMSLevel))
      return ExitWithResult(sResultPath, " Error - invalid ctx_iMSLevel in worker job.\n", 1);

   if (!ReadDoubleListFile(sMassOffsetsFile, ctx.vMassOffsets, sErrorMsg))
      return ExitWithResult(sResultPath, sErrorMsg, 1);

   set<int> setExplicitScans;
   if (bUseExplicitScans)
   {
      if (sExplicitScansFile.empty())
         return ExitWithResult(sResultPath, " Error - explicit scan filter enabled but explicit_scans_file is empty.\n", 1);

      if (!ParseScanIntegersFromFile(sExplicitScansFile,
                                     setExplicitScans,
                                     sErrorMsg,
                                     "prefilter_worker_explicit_scans"))
      {
         return ExitWithResult(sResultPath, sErrorMsg, 1);
      }
   }

   vector<double> vNovelMasses;
   if (bUseNovelMassFilter)
   {
      bool bLoadedFromInternalDetailed = false;
      if (!sInternalNovelPeptideFile.empty())
      {
         vector<NovelPeptideRecord> vInternalRecords;
         vector<double> vPrecomputedMasses;
         bool bHasDetailedMzColumns = false;
         if (!ParseInternalNovelPeptideFile(sInternalNovelPeptideFile,
                                            vInternalRecords,
                                            sErrorMsg,
                                            &vPrecomputedMasses,
                                            &bHasDetailedMzColumns))
         {
            return ExitWithResult(sResultPath, sErrorMsg, 1);
         }

         if (bHasDetailedMzColumns)
         {
            if (vPrecomputedMasses.empty())
               return ExitWithResult(sResultPath, " Error - detailed internal novel peptide file did not contain any mz rows.\n", 1);
            vNovelMasses.swap(vPrecomputedMasses);
            bLoadedFromInternalDetailed = true;
         }
      }

      if (!bLoadedFromInternalDetailed)
      {
         if (sNovelMassesFile.empty())
            return ExitWithResult(sResultPath, " Error - novel mass filter enabled but novel_masses_file is empty.\n", 1);

         if (!ReadDoubleListFile(sNovelMassesFile, vNovelMasses, sErrorMsg))
            return ExitWithResult(sResultPath, sErrorMsg, 1);
      }
   }

   InputFileInfo inputFile;
   inputFile.iAnalysisType = iAnalysisType;
   inputFile.iFirstScan = iFirstScan;
   inputFile.iLastScan = iLastScan;
   strncpy(inputFile.szFileName, sInputPath.c_str(), SIZE_FILE - 1);
   inputFile.szFileName[SIZE_FILE - 1] = '\0';
   inputFile.szBaseName[0] = '\0';

   string sTempMgfPath;
   int iNumScansKept = 0;
   if (!FilterInputFileToTempMgf(inputFile,
                                 setExplicitScans,
                                 bUseExplicitScans,
                                 vNovelMasses,
                                 bUseNovelMassFilter,
                                 ctx,
                                 sTempMgfPath,
                                 iNumScansKept,
                                 sErrorMsg))
   {
      return ExitWithResult(sResultPath, sErrorMsg, 1);
   }

   if (!WriteWorkerResultFile(sResultPath, true, sTempMgfPath, iNumScansKept, ""))
   {
      fprintf(stderr, " Error - failed writing prefilter worker result file \"%s\".\n", sResultPath.c_str());
      return 1;
   }

   return 0;
}
