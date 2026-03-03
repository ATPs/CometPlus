#include "Common.h"
#include "CometData.h"
#include "CometPlusParams.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <map>

using std::set;
using std::string;
using std::vector;

static bool StartsWith(const string &s, const char *pszPrefix)
{
   size_t iLen = strlen(pszPrefix);
   return s.size() >= iLen && s.compare(0, iLen, pszPrefix) == 0;
}

void TrimWhitespace(char *buf)
{
   int iLen = (int)strlen(buf);
   char* szTrimmed = buf;
   while (iLen > 0 && isspace((unsigned char)szTrimmed[iLen - 1]))
      szTrimmed[--iLen] = 0;
   while (*szTrimmed && isspace((unsigned char)*szTrimmed))
   {
      ++szTrimmed;
      --iLen;
   }
   memmove(buf, szTrimmed, iLen + 1);
}

void StripCommentInPlace(char *buf)
{
   char *pComment = strchr(buf, '#');
   if (pComment != NULL)
      *pComment = '\0';
}

const set<string>& GetDedicatedOverrideKeys()
{
   static const set<string> setKeys =
   {
      "database_name",
      "num_threads",
      "scan_range",
      "spectrum_batch_size"
   };

   return setKeys;
}

static const set<string>& GetFallbackHelpParamKeys()
{
   static const set<string> setKeys =
   {
      "activation_method",
      "add_A_alanine",
      "add_B_user_amino_acid",
      "add_C_cysteine",
      "add_Cterm_peptide",
      "add_Cterm_protein",
      "add_D_aspartic_acid",
      "add_E_glutamic_acid",
      "add_F_phenylalanine",
      "add_G_glycine",
      "add_H_histidine",
      "add_I_isoleucine",
      "add_J_user_amino_acid",
      "add_K_lysine",
      "add_L_leucine",
      "add_M_methionine",
      "add_N_asparagine",
      "add_Nterm_peptide",
      "add_Nterm_protein",
      "add_O_pyrrolysine",
      "add_P_proline",
      "add_Q_glutamine",
      "add_R_arginine",
      "add_S_serine",
      "add_T_threonine",
      "add_U_selenocysteine",
      "add_V_valine",
      "add_W_tryptophan",
      "add_X_user_amino_acid",
      "add_Y_tyrosine",
      "add_Z_user_amino_acid",
      "allowed_missed_cleavage",
      "clear_mz_range",
      "clip_nterm_aa",
      "clip_nterm_methionine",
      "correct_mass",
      "database_name",
      "decoy_prefix",
      "decoy_search",
      "digest_mass_range",
      "equal_I_and_L",
      "explicit_deltacn",
      "export_additional_pepxml_scores",
      "fragindex_max_fragmentmass",
      "fragindex_min_fragmentmass",
      "fragindex_min_ions_report",
      "fragindex_min_ions_score",
      "fragindex_num_spectrumpeaks",
      "fragindex_skipreadprecursors",
      "fragment_bin_offset",
      "fragment_bin_tol",
      "isotope_error",
      "mango_search",
      "mass_offsets",
      "mass_type_fragment",
      "mass_type_parent",
      "max_duplicate_proteins",
      "max_fragment_charge",
      "max_iterations",
      "max_precursor_charge",
      "max_variable_mods_in_peptide",
      "min_precursor_charge",
      "minimum_intensity",
      "minimum_peaks",
      "minimum_xcorr",
      "ms1_mass_range",
      "ms_level",
      "nucleotide_reading_frame",
      "num_enzyme_termini",
      "num_output_lines",
      "num_results",
      "num_threads",
      "old_mods_encoding",
      "output_mzidentmlfile",
      "output_pepxmlfile",
      "output_percolatorfile",
      "output_sqtfile",
      "output_sqtstream",
      "output_suffix",
      "output_txtfile",
      "override_charge",
      "peff_format",
      "peff_obo",
      "peff_verbose_output",
      "peptide_length_range",
      "peptide_mass_tolerance",
      "peptide_mass_tolerance_lower",
      "peptide_mass_tolerance_upper",
      "peptide_mass_units",
      "percentage_base_peak",
      "pinfile_protein_delimiter",
      "precursor_NL_ions",
      "precursor_charge",
      "precursor_tolerance_type",
      "print_ascorepro_score",
      "print_expect_score",
      "protein_modslist_file",
      "remove_precursor_peak",
      "remove_precursor_tolerance",
      "require_variable_mod",
      "resolve_fullpaths",
      "sample_enzyme_number",
      "scale_fragmentNL",
      "scan_range",
      "search_enzyme2_number",
      "search_enzyme_number",
      "set_A_alanine",
      "set_B_user_amino_acid",
      "set_C_cysteine",
      "set_D_aspartic_acid",
      "set_E_glutamic_acid",
      "set_F_phenylalanine",
      "set_G_glycine",
      "set_H_histidine",
      "set_I_isoleucine",
      "set_J_user_amino_acid",
      "set_K_lysine",
      "set_L_leucine",
      "set_M_methionine",
      "set_N_asparagine",
      "set_O_pyrrolysine",
      "set_P_proline",
      "set_Q_glutamine",
      "set_R_arginine",
      "set_S_serine",
      "set_T_threonine",
      "set_U_selenocysteine",
      "set_V_valine",
      "set_W_tryptophan",
      "set_X_user_amino_acid",
      "set_Y_tyrosine",
      "set_Z_user_amino_acid",
      "speclib_ms_level",
      "spectral_library_ms_level",
      "spectral_library_name",
      "spectrum_batch_size",
      "text_file_extension",
      "theoretical_fragment_ions",
      "use_A_ions",
      "use_B_ions",
      "use_C_ions",
      "use_NL_ions",
      "use_X_ions",
      "use_Y_ions",
      "use_Z1_ions",
      "use_Z_ions",
      "variable_mod01",
      "variable_mod02",
      "variable_mod03",
      "variable_mod04",
      "variable_mod05",
      "variable_mod06",
      "variable_mod07",
      "variable_mod08",
      "variable_mod09",
      "variable_mod10",
      "variable_mod11",
      "variable_mod12",
      "variable_mod13",
      "variable_mod14",
      "variable_mod15",
      "xcorr_processing_offset"
   };

   return setKeys;
}

static bool IsToggleLikeKey(const string &sKey)
{
   if (StartsWith(sKey, "use_"))
      return true;

   static const set<string> setToggleKeys =
   {
      "clip_nterm_methionine",
      "decoy_search",
      "equal_I_and_L",
      "explicit_deltacn",
      "fragindex_skipreadprecursors",
      "mass_type_fragment",
      "mass_type_parent",
      "old_mods_encoding",
      "output_mzidentmlfile",
      "output_pepxmlfile",
      "output_percolatorfile",
      "output_sqtfile",
      "output_sqtstream",
      "output_txtfile",
      "peff_format",
      "peff_verbose_output",
      "print_expect_score",
      "require_variable_mod",
      "resolve_fullpaths",
      "scale_fragmentNL",
      "theoretical_fragment_ions"
   };

   return setToggleKeys.find(sKey) != setToggleKeys.end();
}

static string BuildFallbackExample(const string &sKey)
{
   if (StartsWith(sKey, "variable_mod"))
      return "15.994915 M 0 3 -1 0 0 0.0";

   if (sKey == "activation_method")
      return "ALL";

   if (sKey == "spectral_library_name" || sKey == "database_name" || sKey == "peff_obo"
         || sKey == "protein_modslist_file")
      return "/path/to/file";

   if (sKey == "output_suffix")
      return ".comet";

   if (sKey == "text_file_extension")
      return "txt";

   if (sKey == "decoy_prefix")
      return "DECOY_";

   if (sKey == "mass_offsets")
      return "0.0";

   if (sKey == "pinfile_protein_delimiter")
      return "tab";

   if (sKey == "activation_method")
      return "ALL";

   if (sKey == "digest_mass_range" || sKey == "scan_range" || sKey == "precursor_charge"
         || sKey == "peptide_length_range" || sKey == "clear_mz_range" || sKey == "ms1_mass_range")
      return "0 0";

   if (sKey == "mass_type_fragment" || sKey == "mass_type_parent")
      return "1";

   if (sKey == "peptide_mass_units")
      return "2";

   if (IsToggleLikeKey(sKey))
      return "0";

   if (StartsWith(sKey, "add_") || StartsWith(sKey, "set_")
         || sKey.find("tol") != string::npos || sKey.find("mass") != string::npos)
      return "0.0";

   if (sKey.find("num_") != string::npos || sKey.find("charge") != string::npos
         || sKey.find("enzyme") != string::npos)
      return "1";

   return "<value>";
}

static string BuildFallbackComment(const string &sKey)
{
   if (StartsWith(sKey, "variable_mod"))
      return "format: <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required> <neutral_loss>";

   if (sKey == "activation_method")
      return "allowed values: ALL, CID, ECD, ETD, ETD+SA, PQD, HCD, IRMPD, SID";

   if (sKey == "decoy_search")
      return "0=no, 1=internal decoy concatenated, 2=internal decoy separate";

   if (sKey == "num_threads")
      return "0=auto detect threads; otherwise explicit thread count";

   if (sKey == "scan_range" || sKey == "digest_mass_range" || sKey == "peptide_length_range"
         || sKey == "clear_mz_range" || sKey == "precursor_charge")
      return "two values: lower upper";

   if (sKey == "spectral_library_ms_level")
      return "alias for internal key speclib_ms_level";

   if (IsToggleLikeKey(sKey))
      return "0/1 flag unless otherwise documented";

   if (sKey.find("file") != string::npos || sKey.find("name") != string::npos)
      return "path/string value";

   return "use params-file value format";
}

const vector<ParamHelpEntry>& GetFallbackHelpEntries()
{
   static const vector<ParamHelpEntry> vEntries = []() -> vector<ParamHelpEntry>
   {
      vector<ParamHelpEntry> vTmp;
      const set<string> &setKeys = GetFallbackHelpParamKeys();
      vTmp.reserve(setKeys.size());

      for (set<string>::const_iterator it = setKeys.begin(); it != setKeys.end(); ++it)
      {
         ParamHelpEntry entry;
         entry.sName = *it;
         entry.sValueExample = BuildFallbackExample(*it);
         entry.sComment = BuildFallbackComment(*it);
         vTmp.push_back(entry);
      }

      return vTmp;
   }();

   return vEntries;
}

bool CollectParamsHelpEntries(const char *pszParamsFile,
                              vector<ParamHelpEntry> &vEntries,
                              string &sErrorMsg)
{
   vEntries.clear();

   FILE *fp = fopen(pszParamsFile, "r");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot open parameter file \"" + string(pszParamsFile) + "\".\n";
      return false;
   }

   char szParamBuf[SIZE_BUF];
   std::map<string, size_t> mapIndexByKey;

   while (!feof(fp))
   {
      if (fgets(szParamBuf, SIZE_BUF, fp) == NULL)
         continue;

      char *pLine = szParamBuf;
      while (*pLine && isspace((unsigned char)*pLine))
         ++pLine;

      if (*pLine == '\0' || *pLine == '#')
         continue;

      if (!strncmp(pLine, "[COMET_ENZYME_INFO]", 19))
         break;

      char *pEquals = strchr(pLine, '=');
      if (pEquals == NULL)
         continue;

      *pEquals = '\0';
      TrimWhitespace(pLine);
      if (*pLine == '\0')
         continue;

      char *pValue = pEquals + 1;
      char *pComment = strchr(pValue, '#');
      if (pComment != NULL)
      {
         *pComment = '\0';
         ++pComment;
         TrimWhitespace(pComment);
      }
      TrimWhitespace(pValue);

      ParamHelpEntry entry;
      entry.sName = pLine;
      entry.sValueExample = pValue;
      entry.sComment = pComment == NULL ? "" : pComment;

      std::map<string, size_t>::iterator it = mapIndexByKey.find(entry.sName);
      if (it == mapIndexByKey.end())
      {
         mapIndexByKey.insert(std::make_pair(entry.sName, vEntries.size()));
         vEntries.push_back(entry);
      }
      else
      {
         vEntries.at(it->second) = entry;
      }
   }

   fclose(fp);
   return true;
}

bool CollectParamsFileKeys(const char *pszParamsFile,
                           set<string> &setParamKeys,
                           string &sErrorMsg)
{
   vector<ParamHelpEntry> vEntries;
   if (!CollectParamsHelpEntries(pszParamsFile, vEntries, sErrorMsg))
      return false;

   setParamKeys.clear();
   for (size_t i = 0; i < vEntries.size(); ++i)
      setParamKeys.insert(vEntries.at(i).sName);

   return true;
}
