#include "ArgParser.h"
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

bool str_compare(const std::string& str1, const std::string& str2) {
	if (str1.size() != str2.size()) { return false; }
	for (std::string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2) {
		if (std::tolower(*c1) != std::tolower(*c2)) { return false;	}
	}
	return true;
}