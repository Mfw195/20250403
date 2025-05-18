function [R_d,R_u,T_d,T_u] = ref_single (Vp_0,Vs_0,rho_0,Vp_1,Vs_1,rho_1,seta)

mu_0 = Vs_0.^2.*rho_0;
mu_1 = Vs_1.^2.*rho_1;
c = 2*(mu_0 - mu_1);

u = 1./Vp_0.*sin(seta);
a0 = sqrt(Vp_0.^(-2) - u.^2);
a1 = sqrt(Vp_1.^(-2) - u.^2);
b0 = sqrt(Vs_0.^(-2) - u.^2);
b1 = sqrt(Vs_1.^(-2) - u.^2);

D1_d = (c.*u.^2 - rho_0 + rho_1).^2.*u.^2 + (c.*u.^2 - rho_0).^2.*a1.*b1 + rho_0.*rho_1.*a1.*b0;
D2_d = c.^2.*u.^2.*a0.*a1.*b0.*b1 + (c.*u.^2 + rho_1).^2.*a0.*b0 + rho_0.*rho_1.*a0.*b1;         %   (6.15)

D1_u = (c.*u.^2 - rho_0 + rho_1).^2.*u.^2 + (c.*u.^2 + rho_1).^2.*a0.*b0 + rho_0.*rho_1.*a0.*b1;
D2_u = c.^2.*u.^2.*a0.*a1.*b0.*b1 + (c.*u.^2 - rho_0).^2.*a1.*b1 + rho_0.*rho_1.*a1.*b0;         %   (6.22)

R_d_pp = (D2_d - D1_d)./(D2_d + D1_d);
R_d_ps = -2*u.*a0./(D2_d + D1_d).*((c.*u.^2 - rho_0 + rho_1).*(c.*u.^2 + rho_1) + c.*(c.*u.^2 - rho_0).*a1.*b1);
T_d_pp = 2*rho_0.*a0./(D2_d + D1_d).*((c.*u.^2 + rho_1).*b0 - (c.*u.^2 - rho_0).*b1);
T_d_ps = -2*rho_0.*u.*a0./(D2_d + D1_d).*(c.*u.^2 - rho_0 + rho_1 + c.*a1.*b0);               %   (6.14)

R_d_sp = 2*u.*b0./(D2_d + D1_d).*((c.*u.^2 - rho_0 + rho_1).*(c.*u.^2 + rho_1) + c.*(c.*u.^2 - rho_0).*a1.*b1);
R_d_ss = (D2_d - D1_d - 2*rho_0.*rho_1.*(a0.*b1 - a1.*b0))./(D1_d + D2_d);
T_d_sp = 2*rho_0.*u.*b0./(D2_d + D1_d).*(c.*u.^2 - rho_0 + rho_1 + c.*a0.*b1);
T_d_ss = 2*rho_0.*b0./(D2_d + D1_d).*((c.*u.^2 + rho_1).*a0 - (c.*u.^2 - rho_0).*a1);         %   (6.18)

R_u_pp = (D2_u - D1_u)./(D2_u + D1_u);
R_u_ps = 2*u.*a1./(D2_u + D1_u).*((c.*u.^2 - rho_0 + rho_1).*(c.*u.^2 - rho_0) + c.*(c.*u.^2 + rho_1).*a0.*b0);
T_u_pp = 2*rho_1.*a1./(D2_u + D1_u).*((c.*u.^2 + rho_1).*b0 - (c.*u.^2 - rho_0).*b1);
T_u_ps = -2*rho_1.*u.*a1./(D2_u + D1_u).*(c.*u.^2 - rho_0 + rho_1 + c.*a0.*b1);               %   (6.21)

R_u_sp = -2*u.*b1./(D2_u + D1_u).*((c.*u.^2 - rho_0 + rho_1).*(c.*u.^2 - rho_0) + c.*(c.*u.^2 + rho_1).*a0.*b0);
R_u_ss = (D2_u - D1_u - 2*rho_0.*rho_1.*(a1.*b0 - a0.*b1))./(D1_u + D2_u);
T_u_sp = 2*rho_1.*u.*b1./(D2_u + D1_u).*(c.*u.^2 - rho_0 + rho_1 + c.*a1.*b0);
T_u_ss = 2*rho_1.*b1./(D2_u + D1_u).*((c.*u.^2 + rho_1).*a0 - (c.*u.^2 - rho_0).*a1);         %   (6.25)

R_d = [ R_d_pp , R_d_sp ; R_d_ps , R_d_ss ];
R_u = [ R_u_pp , R_u_sp ; R_u_ps , R_u_ss ];
T_d = [ T_d_pp , T_d_sp ; T_d_ps , T_d_ss ];
T_u = [ T_u_pp , T_u_sp ; T_u_ps , T_u_ss ];                                  %   （6.35）








