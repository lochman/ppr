#include <vector>
#include "approx\src\CommonApprox.h"

class ArgParser {
	std::vector <std::string> args;
public:
	ArgParser(int argc, char **argv);
	std::string get_option(const std::string &option);
	bool check_option(const std::string &option);
	std::string get_segment();
};

bool str_compare(const std::string &str1, const std::string &str2);
HRESULT parse_method(const std::string &method, IGlucoseLevels *levels, CCommonApprox **approx);