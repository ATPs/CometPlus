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

#ifndef _COMETPLUS_NOVELMODEUTILS_H_
#define _COMETPLUS_NOVELMODEUTILS_H_

#include "Common.h"
#include "CometData.h"
#include "CometInterfaces.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <utility>

struct NovelModeOptions
{
   string sNovelProteinPath;
   string sNovelPeptidePath;
   string sScanFilePath;
   string sScanNumbers;
   string sOutputFolder;
   string sOutputInternalNovelPeptidePath;
   string sInternalNovelPeptidePath;
   bool bStopAfterSavingNovelPeptide = false;

   bool HasNovelInputOptions() const
   {
      return !sNovelProteinPath.empty() || !sNovelPeptidePath.empty();
   }

   bool HasInternalNovelInput() const
   {
      return !sInternalNovelPeptidePath.empty();
   }

   bool HasNovelMode() const
   {
      return HasNovelInputOptions() || HasInternalNovelInput();
   }

   bool HasExplicitScanFilter() const
   {
      return !sScanFilePath.empty() || !sScanNumbers.empty();
   }
};

struct PeptideMassEntry
{
   string sPeptide;
   double dMass;
   bool bDecoy;
   vector<string> vProteinIds;
   vector<unsigned char> vVarModSites;
};

struct NovelPeptideRecord
{
   string sPeptide;
   string sPeptideId;
   vector<string> vProteinIds;
};

struct NovelMassFilterContext
{
   double dPepMassLow;
   double dPepMassHigh;
   double dTolLower;
   double dTolUpper;
   int iTolUnits;
   int iTolType;
   int iIsotopeError;
   vector<double> vMassOffsets;

   int iStartCharge;
   int iEndCharge;
   int iOverrideCharge;
   int iMinPrecursorCharge;
   int iMaxPrecursorCharge;
   bool bCorrectMass;
   int iMSLevel;
};

bool IsIdxDatabasePath(const string& sPath);
bool BuildMergedFasta(const vector<string>& vDatabases,
                      string& sMergedDatabase,
                      string& sErrorMsg);
bool ConfigureDatabaseInputs(const vector<string>& vDatabases,
                             CometInterfaces::ICometSearchManager* pSearchMgr,
                             string& sMergedDatabasePath,
                             string& sErrorMsg);
bool BuildTemporaryParamsFile(const string& sBaseParamsPath,
                              const map<string, string>& mOverrides,
                              string& sTmpParamsPath,
                              string& sErrorMsg);
bool ParsePeptideIdxEntries(const string& sIdxPath,
                            vector<PeptideMassEntry>& vEntries,
                            bool bLoadProteinIds,
                            string& sErrorMsg);
bool ParseFragmentIdxEntries(const string& sIdxPath,
                             vector<PeptideMassEntry>& vEntries,
                             string& sErrorMsg);
bool ParseScanIntegersFromString(const string& sInput,
                                 set<int>& setScans,
                                 string& sErrorMsg,
                                 const string& sSourceLabel);
bool ParseScanIntegersFromFile(const string& sPath,
                               set<int>& setScans,
                               string& sErrorMsg,
                               const string& sSourceLabel);
bool ParseNovelPeptideFile(const string& sPath,
                           vector<string>& vPeptides,
                           string& sErrorMsg);
bool FindNoCutEnzymeNumber(const string& sParamsFilePath,
                           int& iNoCutEnzyme,
                           string& sErrorMsg);
bool RunCometForIndexGeneration(const string& sExecutablePath,
                                const string& sParamsPath,
                                const string& sDatabasePath,
                                bool bFragmentIndex,
                                int iThreadOverride,
                                string& sErrorMsg);
void SetCometPlusTempDirectory(const string& sTempDir);
bool EnsureDirectoryExistsRecursive(const string& sDirPath,
                                    string& sErrorMsg);
bool CreateTempPath(const string& sPrefix,
                    const string& sSuffix,
                    string& sOutPath,
                    string& sErrorMsg);
bool WritePeptidesToFasta(const vector<string>& vPeptides,
                          const string& sPrefix,
                          string& sOutPath,
                          string& sErrorMsg);
bool WriteNovelRecordsToFasta(const vector<NovelPeptideRecord>& vRecords,
                              const string& sPrefix,
                              string& sOutPath,
                              string& sErrorMsg);
bool ParseInternalNovelPeptideFile(const string& sPath,
                                   vector<NovelPeptideRecord>& vRecords,
                                   string& sErrorMsg,
                                   vector<double>* pvPrecomputedMasses = NULL,
                                   bool* pbHasDetailedMzColumns = NULL);
bool BuildNovelMassFilterContext(CometInterfaces::ICometSearchManager* pSearchMgr,
                                 NovelMassFilterContext& ctx,
                                 string& sErrorMsg);
bool FilterInputFileToTempMgf(const InputFileInfo& inputFile,
                              const set<int>& setExplicitScans,
                              bool bUseExplicitScans,
                              const vector<double>& vNovelMasses,
                              bool bUseNovelMassFilter,
                              const NovelMassFilterContext& ctx,
                              string& sOutTempMgfPath,
                              int& iNumScansKept,
                              string& sErrorMsg);
bool BuildNovelMergedSourceLabel(size_t iOneBasedIndex,
                                 const string& sInputStem,
                                 string& sOutLabel,
                                 string& sErrorMsg);
bool RewriteMgfTitleWithSourceLabel(const string& sTitleLine,
                                    const string& sSourceLabel,
                                    string& sOutTitleLine);
bool MergeFilteredMgfFilesWithSourceTag(
      const vector<std::pair<string, string>>& vFilteredMgfAndSourceLabel,
      string& sOutMergedMgfPath,
      string& sErrorMsg);
string NormalizePeptideForCompare(const string& sPeptide, bool bTreatSameIL);
string NormalizePeptideToken(const string& sToken);
string ComputeInputBaseName(const string& sInputPath);

#endif
