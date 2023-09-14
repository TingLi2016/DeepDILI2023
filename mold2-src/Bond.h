#pragma once

class Bond
{
public:
	Bond(void);
	~Bond(void);
	static const char SINGLE_BOND = 1;
	static const char DOUBLE_BOND = 2;
	static const char TRIPLE_BOND = 3;
	static const char AROMATIC_BOND = 4;
	static const char CHAIN_BOND = 2;

private:
	short m_first_atom_number;
	short m_second_atom_number;
	char m_bond_type;
	char m_bond_stereo;
	char m_bond_topology;
	char m_reacting_center_status;
public:
	void setFirstAtomNumber(short nAtomID);
	short getFirstAtomNumber();
	void setSecondAtomNumber(short nAtomID);
	short getSecondAtomNumber();
	void setBondType(char bondType);
	char getBondType();
	void setBondStereo(char bondStereo);
	char getBondStereo();
	void setBondTopology(char bondTopology);
	char getBondTopology();
	void setReactingCenterStatus(char reactingCenterStatus);
	char getReactingCenterStatus();
};
