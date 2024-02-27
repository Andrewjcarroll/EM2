const double amp1 = em3::EM3_ID_AMP1;
const double lambda1 = em3::EM3_ID_LAMBDA1;

double B0, B1, B2, E0, E1, E2;
double rho_e, J0, J1, J2;
double r, Ephi_over_r, tmp_Ephiup;

r = sqrt(x * x + y * y + z * z);
tmp_Ephiup = -8.0 * amp1 * lambda1 * lambda1 * exp(-lambda1 * r * r);
E0 = -y * tmp_Ephiup;
E1 = x * tmp_Ephiup;
E2 = 0.0;

B0 = 0.0;
B1 = 0.0;
B2 = 0.0;

J0 = 0.0;
J1 = 0.0;
J2 = 0.0;

rho_e = 0.0;

var[VAR::U_E0] = E0;
var[VAR::U_E1] = E1;
var[VAR::U_E2] = E2;
var[VAR::U_B0] = B0;
var[VAR::U_B1] = B1;
var[VAR::U_B2] = B2;

