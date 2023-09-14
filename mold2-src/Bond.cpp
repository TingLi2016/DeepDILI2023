#include "stdafx.h"
#include "Bond.h"

Bond::Bond(void)
: m_first_atom_number(0), m_second_atom_number(0), m_bond_type(1), m_bond_stereo(0),
m_bond_topology(0), m_reacting_center_status(0)
{
}

Bond::~Bond(void)
{
}

void Bond::setFirstAtomNumber(short nAtomID)
{
	m_first_atom_number = nAtomID;
}

short Bond::getFirstAtomNumber()
{
	return m_first_atom_number;
}


void Bond::setSecondAtomNumber(short nAtomID)
{
	m_second_atom_number = nAtomID;
}

short Bond::getSecondAtomNumber()
{
	return m_second_atom_number;
}

void Bond::setBondType(char bondType)
{
	m_bond_type = bondType;
}

char Bond::getBondType()
{
	return m_bond_type;
}

void Bond::setBondStereo(char bondStereo)
{
	m_bond_stereo = bondStereo;
}

char Bond::getBondStereo()
{
	return m_bond_stereo;
}

void Bond::setBondTopology(char bondTopology)
{
	m_bond_topology = bondTopology;
}

char Bond::getBondTopology()
{
	return m_bond_topology;
}


void Bond::setReactingCenterStatus(char reactingCenterStatus)
{
	m_reacting_center_status = reactingCenterStatus;
}

char Bond::getReactingCenterStatus()
{
	return m_reacting_center_status;
}
