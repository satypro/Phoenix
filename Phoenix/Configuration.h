#pragma once
#include<map>

class Configuration
{
public:
	Configuration();
	~Configuration();
	void SetValue(std::string key, std::string value);
	std::string GetValue(std::string key);
	std::map<std::string, std::string>& GetConfig();
private:
	std::map<std::string, std::string> _config;
};

