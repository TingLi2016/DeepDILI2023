#pragma once

class Constants
{
public:
	static const short ERROR_SIGNAL;
	static const short MAX_CONNECTIONS;
	static const double ROUNDOFF;
	static const float LOG20;

	static const float CovalentRadius1[];
	static const float CovalentRadius[];
	static const double val_ele_pri_num[];
	static const float polarizability_value[];
	static const signed char defined_charge[];
	static const char atom_valences[][5];
	static const float atom_weight[];
	static const float vdw_radius[];
	static const float vdw_volume[];
	static const float elect_negativities[][3];

	Constants(void);
	~Constants(void);
};
