t1 = gd[0][0];
t2 = gd[1][1];
t4 = gd[2][2];
t6 = gd[1][2];
t7 = t6 * t6;
t9 = gd[0][1];
t10 = t9 * t9;
t12 = gd[0][2];
t16 = t12 * t12;
detgd = t1 * t2 * t4 - t1 * t7 - t10 * t4 + 2.0 * t9 * t12 * t6 - t16 * t2;
idetgd = 0.1E1 / detgd;
gu[0][0] = idetgd * (gd[1][1] * gd[2][2] - gd[1][2] * gd[1][2]);
gu[0][1] = idetgd * (-gd[0][1] * gd[2][2] + gd[0][2] * gd[1][2]);
gu[0][2] = idetgd * (gd[0][1] * gd[1][2] - gd[0][2] * gd[1][1]);
gu[1][0] = idetgd * (-gd[0][1] * gd[2][2] + gd[0][2] * gd[1][2]);
gu[1][1] = idetgd * (gd[0][0] * gd[2][2] - gd[0][2] * gd[0][2]);
gu[1][2] = idetgd * (-gd[0][0] * gd[1][2] + gd[0][1] * gd[0][2]);
gu[2][0] = idetgd * (gd[0][1] * gd[1][2] - gd[0][2] * gd[1][1]);
gu[2][1] = idetgd * (-gd[0][0] * gd[1][2] + gd[0][1] * gd[0][2]);
gu[2][2] = idetgd * (gd[0][0] * gd[1][1] - gd[0][1] * gd[0][1]);

chi = pow(idetgd, 0.3333333333333333);

gtd[0][0] = chi * gd[0][0];
gtd[0][1] = chi * gd[0][1];
gtd[0][2] = chi * gd[0][2];
gtd[1][0] = chi * gd[0][1];
gtd[1][1] = chi * gd[1][1];
gtd[1][2] = chi * gd[1][2];
gtd[2][0] = chi * gd[0][2];
gtd[2][1] = chi * gd[1][2];
gtd[2][2] = chi * gd[2][2];

trK = gu[0][0] * Kd[0][0] + 2.0 * gu[0][1] * Kd[0][1] +
      2.0 * gu[0][2] * Kd[0][2] + gu[1][1] * Kd[1][1] +
      2.0 * gu[1][2] * Kd[1][2] + gu[2][2] * Kd[2][2];

Atd[0][0] = chi * (Kd[0][0] - 0.3333333333333333 * gd[0][0] * trK);
Atd[0][1] = chi * (Kd[0][1] - 0.3333333333333333 * gd[0][1] * trK);
Atd[0][2] = chi * (Kd[0][2] - 0.3333333333333333 * gd[0][2] * trK);
Atd[1][0] = Atd[0][1];
Atd[1][1] = chi * (Kd[1][1] - 0.3333333333333333 * gd[1][1] * trK);
Atd[1][2] = chi * (Kd[1][2] - 0.3333333333333333 * gd[1][2] * trK);
Atd[2][0] = Atd[0][2];
Atd[2][1] = Atd[1][2];
Atd[2][2] = chi * (Kd[2][2] - 0.3333333333333333 * gd[2][2] * trK);
