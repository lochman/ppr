#include <vector>

class ArgParser {
	std::vector <std::string> args;
public:
	ArgParser(int argc, char **argv);
	std::string get_option(const std::string &option);
	bool check_option(const std::string &option);
};

bool str_compare(const std::string& str1, const std::string& str2);