#pragma once

class mold2Exception :
	public std::exception
{
private:
	std::string m_msg;
  std::string m_molid;
public:
	const char* what() const throw(){
		return m_msg.c_str();
	}
	const char* molid() const throw(){
		return m_molid.c_str();
	}

	mold2Exception(std::string msg): m_msg(msg), m_molid("") { }
	mold2Exception(std::string msg, std::string molId): m_msg(msg), m_molid(molId) { }

	~mold2Exception(void) throw()
	{
	}
};
