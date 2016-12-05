#include "ArgParser.h"
#include "approx\src\CubicApprox.h"
#include "approx\src\AkimaApprox.h"
#include "approx\src\CatmullRomApprox.h"
#include <cctype>

ArgParser::ArgParser(int argc, char **argv) {
	for (int i = 1; i < argc; ++i) { args.push_back(std::string(argv[i])); }
}

std::string ArgParser::get_option(const std::string &option) {
	std::vector<std::string>::const_iterator itr = std::find(args.begin(), args.end(), option);
	if (itr != args.end() && ++itr != args.end()) { return *itr; }
	return "";
}

bool ArgParser::check_option(const std::string &option) {
	return std::find(args.begin(), args.end(), option) != args.end();
}

bool str_compare(const std::string &str1, const std::string &str2) {
	if (str1.size() != str2.size()) { return false; }
	for (std::string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2) {
		if (std::tolower(*c1) != std::tolower(*c2)) { return false;	}
	}
	return true;
}

HRESULT parse_method(const std::string &method, IGlucoseLevels *levels, CCommonApprox **approx) {
	if (str_compare(method, "akima") || str_compare(method, "a")) {
		*approx = new AkimaApprox(levels);
	} else if (str_compare(method, "catmull") || str_compare(method, "cr")) {
		*approx = new CatmullRomApprox(levels);
	} else if (str_compare(method, "cubic") || str_compare(method, "c")) {
		*approx = new CubicApprox(levels);
	} else {
		*approx = new AkimaApprox(levels);
	}
	return S_OK;
}