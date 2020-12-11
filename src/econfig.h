#ifndef _ECONFIG_H_
#define _ECONFIG_H_

typedef struct tEgulpConfigure {
  double dh1sz;
  double dh2sz;
  double dh3sz;
  double vdw_factor_i;
  double vdw_factor_f;
  int use_vdw_factor;
  double offset;
  int build_grid;
  int build_grid_from_scratch;
  int save_grid;
  char input_grid_file[128];
  char output_grid_file[128];
  char output_pot_file[128];
  int calculate_pot_diff;
  int calculate_pot;
  int skip_everything;
  int point_charges_present;
  int include_pceq;
  int imethod;
} tConfigure;

void InitEconfigDefaults(tConfigure *confptr);
void InitEconfig(tConfigure *confptr, char *fname);

void DestroyEconfig(tConfigure conf);

#endif
