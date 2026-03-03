#ifndef COMETPLUS_PARAMS_H
#define COMETPLUS_PARAMS_H

#include <set>
#include <string>
#include <vector>

namespace CometInterfaces
{
   class ICometSearchManager;
}

struct CmdParamOverride
{
   std::string sName;
   std::string sValue;
};

struct ParamHelpEntry
{
   std::string sName;
   std::string sValueExample;
   std::string sComment;
};

void TrimWhitespace(char *buf);
void StripCommentInPlace(char *buf);

const std::set<std::string>& GetDedicatedOverrideKeys();
const std::vector<ParamHelpEntry>& GetFallbackHelpEntries();

bool CollectParamsFileKeys(const char *pszParamsFile,
                           std::set<std::string> &setParamKeys,
                           std::string &sErrorMsg);
bool CollectParamsHelpEntries(const char *pszParamsFile,
                              std::vector<ParamHelpEntry> &vEntries,
                              std::string &sErrorMsg);

void LoadParameters(char *pszParamsFile,
                    CometInterfaces::ICometSearchManager *pSearchMgr,
                    const std::vector<CmdParamOverride> &vCliParamOverrides);

#endif
