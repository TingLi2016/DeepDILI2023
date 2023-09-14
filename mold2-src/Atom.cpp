#include "stdafx.h"
#include "Atom.h"

Atom::Atom(void)
: m_x(0),m_y(0),m_z(0),m_atom_symbol(0),m_mass_difference(0), m_charge(0),m_stereo_Parity(0),
m_hydrogen_count_plus1(0), m_stereo_care_box(0), m_valence(0), m_H0_designator(0), 
m_reaction_component_type(0), m_reaction_component_number(0), 
m_atom_atom_mapping_number(0), m_inversion_retention_Flag(0), m_exact_change_Flag(0){

}

Atom::~Atom(void)
{
}

void Atom::setAtomSymbol(char symbol)
{
	m_atom_symbol = symbol;
}

char Atom::getAtomSymbol()
{
	return m_atom_symbol;
}

void Atom::setX(float x)
{
	m_x = x;
}

void Atom::setY(float y)
{
	m_y = y;
}

void Atom::setZ(float z)
{
	m_z = z;
}

void Atom::setCharge(char charge)
{
	m_charge = charge;
}

void Atom::setAtomMappingNumber(char atomMappingNumber)
{
	this->m_atom_atom_mapping_number = atomMappingNumber;
}

char Atom::getAtomMappingNumber(void)
{
	return this->m_atom_atom_mapping_number;
}

void Atom::setValence(char valence)
{
	m_valence = valence;
}

void Atom::setExactChangeFlag(char exactChangeFlag)
{
	m_exact_change_Flag = exactChangeFlag;
}

void Atom::setH0designator(char H0designator)
{
	m_H0_designator = H0designator;
}

void Atom::setHydrogenCountPlus1(char hydrogenCountPlus1)
{
	m_hydrogen_count_plus1 = hydrogenCountPlus1;
}

void Atom::setMassDifference(char massDifference)
{
	m_mass_difference = massDifference;
}

void Atom::setStereoCareBox(char stereoCareBox)
{
	m_stereo_care_box = stereoCareBox;
}

void Atom::setStereoParity(char stereoParity)
{
	m_stereo_Parity = stereoParity;
}

void Atom::setReactionComponentNumber(char reactionComponentNumber)
{
	m_reaction_component_number = reactionComponentNumber;
}

void Atom::setReactionComponentType(char reactionComponentType)
{
	m_reaction_component_type = reactionComponentType;
}

void Atom::setInversionRetentionFlag(char inversionRetensionFlag)
{
	m_inversion_retention_Flag = inversionRetensionFlag;
}

char Atom::getHydrogenCountPlus1()
{
	return m_hydrogen_count_plus1;
}

char Atom::getValence(void)
{
	return m_valence;
}

char Atom::getCharge(void)
{
	return m_charge;
}

float Atom::getX(void)
{
	return m_x;
}

float Atom::getY(void)
{
	return m_y;
}

float Atom::getZ(void)
{
	return m_z;
}

char Atom::getMassDifference(void)
{
	return this->m_mass_difference;
}

char Atom::getAtomStereoParity(void)
{
	return this->m_stereo_Parity;
}

char Atom::getH0_designator(void)
{
	return this->m_H0_designator;
}

char Atom::getReactionComponentType(void)
{
	return this->m_reaction_component_type;
}

char Atom::getReactionComponentNumber(void)
{
	return this->m_reaction_component_number;
}

char Atom::getInversionRetentionFlag(void)
{
	return this->m_inversion_retention_Flag;
}

char Atom::getExactChangeFlag(void)
{
	return this->m_exact_change_Flag;
}

char Atom::getStereoCareBox(void)
{
	return this->m_stereo_care_box;
}
