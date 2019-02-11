#include "Configuration.h"

Configuration::Configuration()
{
}

Configuration::~Configuration()
{
}

void Configuration::SetValue(std::string key, std::string value)
{
	std::map<std::string, std::string>::iterator it;

	it = _config.find(key);

	if (it != _config.end())
	{
		it->second = value;
	}
	else
	{
		_config.insert(std::make_pair(key, value));
	}
}

std::string Configuration::GetValue(std::string key)
{
	std::map<std::string, std::string>::iterator it;

	it = _config.find(key);
	if (it != _config.end())
		return it->second;
	return "";
}

std::map<std::string, std::string>& Configuration::GetConfig()
{
	return _config;
}