// Copyright 2026 Jimmy Eng
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


#ifndef _COMETPLUSMULTIDB_H_
#define _COMETPLUSMULTIDB_H_

#include "Common.h"
#include "CometDataInternal.h"

struct CometPlusIdxContext
{
   string sPath;
   int iIndexDbType;  // 1=fragment index, 2=peptide index

   // Local peptide/protein-set entry index -> global protein-set entry index.
   vector<int> vLocalProteinSetToGlobalSet;

   // Local protein-name offset in this .idx -> global protein id.
   map<comet_fileoffset_t, int> mLocalProteinOffsetToGlobalId;
};

extern bool g_bCometPlusMultiDbMode;
extern bool g_bCometPlusMultiIdxMode;
extern vector<string> g_vCometPlusDatabaseList;
extern vector<CometPlusIdxContext> g_vCometPlusIdxContexts;

// Stable global protein registries used in multi-idx mode.
extern vector<string> g_vCometPlusGlobalProteinNames;
extern vector<vector<int>> g_vCometPlusGlobalProteinSets;

void CometPlusResetMultiDBState();
vector<string> CometPlusSplitDatabaseList(const string& sList);

bool CometPlusProbeIdxType(const string& sPath, int& iIndexDbType, string& sErrorMsg);
bool CometPlusInitializeMultiIdxContexts(const vector<string>& vDbList,
                                         int iExpectedIndexDbType,
                                         string& sErrorMsg);

bool CometPlusReadProteinSetsForContext(FILE* fp,
                                        comet_fileoffset_t clProteinsFilePos,
                                        CometPlusIdxContext& ctx,
                                        string& sErrorMsg);

int CometPlusMergeProteinSetIndices(comet_fileoffset_t lSetA, comet_fileoffset_t lSetB);
bool CometPlusGetProteinSet(comet_fileoffset_t lSetIndex, vector<int>& vOutProteinIds);
bool CometPlusGetFirstProteinIdInSet(comet_fileoffset_t lSetIndex, int& iFirstProteinId);
bool CometPlusGetProteinNameById(int iProteinId, string& sName);

#endif // _COMETPLUSMULTIDB_H_
