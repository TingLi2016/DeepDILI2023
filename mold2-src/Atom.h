#pragma once

class Atom
{
public:
	Atom(void);
	~Atom(void);
private:
	float m_x, m_y, m_z;
	char m_atom_symbol;
	char m_mass_difference;
	char m_charge;
	char m_stereo_Parity;
	char m_hydrogen_count_plus1;
	char m_stereo_care_box;
	char m_valence;
	char m_H0_designator;
	char m_reaction_component_type;
	char m_reaction_component_number;
	char m_atom_atom_mapping_number;
	char m_inversion_retention_Flag;
	char m_exact_change_Flag;
public:
	void setAtomSymbol(char symbol);
	char getAtomSymbol();
	void setX(float x);
	void setY(float y);
	void setZ(float z);
	void setCharge(char charge);
	void setAtomMappingNumber(char atomMappingNumber);
	void setValence(char valence);
	void setExactChangeFlag(char exactChangeFlag);
	void setH0designator(char H0designatro);
	void setHydrogenCountPlus1(char hydrogenCountPlus1);
	void setMassDifference(char massDifference);
	void setStereoCareBox(char stereoCareBox);
	void setStereoParity(char stereoParity);
	void setReactionComponentNumber(char reactionComponentNumber);
	void setReactionComponentType(char reactionComponentType);
	void setInversionRetentionFlag(char inversionRetensionFlag);

	char getHydrogenCountPlus1();
	char getValence(void);
	char getCharge(void);
	float getX(void);
	float getY(void);
	float getZ(void);
	char getMassDifference(void);
	char getAtomStereoParity(void);
	char getH0_designator(void);
	char getReactionComponentType(void);
	char getReactionComponentNumber(void);
	char getInversionRetentionFlag(void);
	char getExactChangeFlag(void);
	char getStereoCareBox(void);
	char getAtomMappingNumber(void);
};
