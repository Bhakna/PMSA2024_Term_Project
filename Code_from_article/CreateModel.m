
% Setup the compartment concentrations

m.Add(model.Concentration( 's', data.V.s ));
m.Add(model.Concentration( 'ri', data.V.ri ));
m.Add(model.Concentration( 'po', data.V.po ));
m.Add(model.Concentration( 'adi', data.V.adi ));
m.Add(model.Concentration( 'h', data.V.h ));
m.Add(model.Concentration( 'g', data.V.g ));
m.Add(model.Concentration( 'art', data.V.art ));
m.Add(model.Concentration( 'vein', data.V.vein ));

m.Add(model.Concentration( 'agut', 1, [], ET_dose_nmol ));

m.Add(model.Concentration( 'clear', 1 ));


% Setup the one dimensional concentration states

A_elf_u = data.N*pi.*((data.D/2).^2 -((data.D/2) - data.h_tb).^2);
A_ep_u  = data.N*pi.*(((data.D/2) + data.h_ep_u).^2-(data.D/2).^2);
A_sub_u = data.N*pi.*(((data.D/2) + data.h_ep_u + data.h_sub_u).^2-((data.D/2) + data.h_ep_u).^2);

A_elf_l = data.N*pi.*((data.D/2).^2 -((data.D/2) - data.h_al).^2);
A_ep_l  = data.N*pi.*(((data.D/2) + data.h_ep_l).^2-(data.D/2).^2);
A_sub_l = data.N*pi.*(((data.D/2) + data.h_ep_l + data.h_sub_l).^2-((data.D/2) + data.h_ep_l).^2);



f_u   = model.Concentration( 'f_u', A_elf_u.AsVolumes(x_upper), x_upper );
ep_u  = model.Concentration( 'ep_u', A_ep_u.AsVolumes(x_upper), x_upper );
sub_u = model.Concentration( 'sub_u', A_sub_u.AsVolumes(x_upper), x_upper );

% (scale alv)

A_area  = data.N*pi.*data.D.* data.L ;
A_alv_tot = A_area.Integral(x_lower);
A_elf_l = A_elf_l + data.h_al   .* data.alv_frac * (data.A_alv_lit - A_alv_tot) ./ data.L;
A_ep_l  = A_ep_l  + data.h_ep_l .* data.alv_frac * (data.A_alv_lit - A_alv_tot) ./ data.L;

f_l   = model.Concentration( 'f_l', A_elf_l.AsVolumes(x_lower), x_lower );
ep_l  = model.Concentration( 'ep_l', A_ep_l.AsVolumes(x_lower), x_lower );
sub_l = model.Concentration( 'sub_l', A_sub_l.AsVolumes(x_lower), x_lower );

% (scale sub)

sub_tot = data.V.lung - f_u.Volume() - ep_u.Volume() - f_l.Volume() - ep_l.Volume();

frac = sub_u.Volume() / sub_l.Volume();
sub_u.ScaleToVolume( frac/(frac+1) * sub_tot );
sub_l.ScaleToVolume( 1/(frac+1) * sub_tot );

% Finally, add the newly created variables
m.Add([f_l, ep_l, sub_l,f_u, ep_u, sub_u]);


% Internal lung flows

D_sub_l = data.D + 2*data.h_ep_l;
D_sub_u = data.D + 2*data.h_ep_u;
  
area_alv = data.alv_frac * (data.A_alv_lit - A_alv_tot) ./ data.L;

q_ep_u = pi*( data.N.*data.D ).*drug_data.P;
q_sub_u = pi*( data.N.*D_sub_u ).*drug_data.P./drug_data.kp_u;

q_ep_l = pi*( data.N.*data.D + area_alv ).*drug_data.P;
q_sub_l = pi*( data.N.*D_sub_l + area_alv ).*drug_data.P./drug_data.kp_u;

dx_l = sub_l.Length();
dx_u = sub_u.Length();

m.Add( model.DualFlow( 'f_l', 'ep_l', ...
        drug_data.fu_f*q_ep_l.OnGrid(x_lower).*dx_l, q_ep_l.OnGrid(x_lower)/drug_data.kp_u.*dx_l ));
    
m.Add( model.DualFlow( 'f_u', 'ep_u', ...
        drug_data.fu_f*q_ep_u.OnGrid(x_upper).*dx_u, q_ep_u.OnGrid(x_upper)/drug_data.kp_u.*dx_u ));

m.Add( model.DualFlow( 'ep_l', 'sub_l', ...
        q_sub_l.OnGrid(x_lower).*dx_l, q_sub_l.OnGrid(x_lower).*dx_l ));
m.Add( model.DualFlow( 'ep_u', 'sub_u', ...
        q_sub_u.OnGrid(x_upper).*dx_u, q_sub_u.OnGrid(x_upper).*dx_u ));

% Connect PBPK to lung

% 2a) Use empiric relationship from Bernard et al. 1996 in lower lung
F = 0.19 + 2.8*exp(-0.51*data.D.OnGrid(x_upper)*1e2);
A_sub_u_F    = F.*A_sub_u.OnGrid(x_upper);

q_br = (sub_u.vol.*F) / sum(sub_u.vol.*F) * data.Q.br;
q_co = sub_l.vol / sub_l.Volume() * data.Q_co;

m.Add( model.Flow( 'vein', 'sub_l', q_co ));
m.Add( model.Flow( 'sub_l', 'art', drug_data.b_p/drug_data.k_p(7)*q_co ));

m.Add( model.Flow( 'art', 'sub_u', q_br ));
m.Add( model.Flow( 'sub_u', 'vein', drug_data.b_p/drug_data.k_p(7)*q_br ));


% Add clearance
m.Add( model.Flow( 'vein', 'clear', drug_data.Cl ));
m.Add( model.Flow( 'agut', 'g', drug_data.F*drug_data.ka ));
m.Add( model.Flow( 'agut', 'clear', (1-drug_data.F)*drug_data.ka ));


% Setup PBPK flow
source = {'s','ri','po','adi','h','g'};
target = {'h','vein','vein','vein','vein','h'};
for k = 1:6
    Q = data.Q.(source{k});
    if k == 5
        Q = Q + data.Q.(source{1}) + data.Q.(source{6});
    end
    m.Add( model.Flow(  source{k}, target{k}, Q*drug_data.b_p / drug_data.k_p(k)));
end

source = 'art';
target = {'s','ri','po','adi','h','g'};
for k = 1:6
    Q = data.Q.(target{k});
    m.Add( model.Flow(  source, target{k}, Q ));
end

% Finally, we add dissolution to our model
m.Add( model.Particles( 'p_l', x_lower, r.^2, drug_data.dens, p_0_l ) );
m.Add( model.Particles( 'p_u', x_upper, r.^2, drug_data.dens, p_0_u ) );

m.Add( model.Dissolution( 'p_l', 'f_l', data.D, drug_data ));
m.Add( model.Dissolution( 'p_u', 'f_u', data.D, drug_data, init_transport_velocity, 'agut' ));
