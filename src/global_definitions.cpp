#include "global_definitions.h"

std::string strtoken(std::string &in, std::string break_symbs)
{
	std::string out;
	while (!in.empty())
	{
		char a = in.front();
		in.erase(in.begin());
		bool break_ = false;
		for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
			if (a == *h) {
				break_ = true;
				break;
			}
		if ((break_) && (out.empty()))
			continue;
		if (break_)
			return out;
		out.push_back(a);
	}
	return out;
}

void ensure_file(std::string fname)
{
	std::string folder = fname;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
#if defined(__WIN32__)
	if (!folder.empty()) {
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
			int code = system(("mkdir \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << GetLastError() << std::endl;
		}
	}
#else
	struct stat st;
	stat(folder.c_str(), &st);
	if (!S_ISDIR(st.st_mode)) {
		int code = system(("mkdir -p \"" + folder + "\"").c_str());
		if (code)
			std::cout << "mkdir -p error: " << code << std::endl;
	}
#endif //_WIN32__
}

char* c_str_cp (const std::string &str)
{
	std::size_t i_end_= str.size();
	char* ret = new char [i_end_+1];
	for (std::size_t i=0; i!=i_end_; ++i) {
		ret[i] = str[i];
	}
	ret[i_end_] = NULL;
	return ret;
}
