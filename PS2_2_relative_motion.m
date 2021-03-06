% known position vectors
r_moon_sc = +140000;
r_sun_earth = -149597898;
r_earth_moon = +384400;
r_sun_jupiter = +778279959;

% known Gm values
Gm_earth = 398600.4415;
Gm_sun = 132712440017.99;
Gm_moon = 4902.8005821478;
Gm_jupiter = 126712767.8578;

fprintf('Acceleration of sc relative to moon:\n');
moon_sc(Gm_earth, Gm_jupiter, Gm_moon, Gm_sun, ...
r_moon_sc, r_sun_earth, r_earth_moon, r_sun_jupiter)

fprintf('Acceleration of sc relative to earth:\n');
earth_sc(Gm_earth, Gm_jupiter, Gm_moon, Gm_sun, ...
r_moon_sc, r_sun_earth, r_earth_moon, r_sun_jupiter)


% spacecraft relative to moon
function moon_sc(Gm_earth, Gm_jupiter, Gm_moon, Gm_sun, ...
    r_moon_sc, r_sun_earth, r_earth_moon, r_sun_jupiter)

    % position vectors relative to spacecraft (sc)
    r_sc_earth = -r_earth_moon - r_moon_sc;
    r_sc_sun = -r_sun_earth - r_earth_moon - r_moon_sc;
    r_sc_jupiter = r_sun_jupiter - r_sun_earth - r_earth_moon - r_moon_sc;
    fprintf('r_sc_earth = %.5e\n', r_sc_earth);
    fprintf('r_sc_sun = %.5e\n', r_sc_sun);
    fprintf('r_sc_jupiter = %.5e\n', r_sc_jupiter);

    % position vectors relative to moon
    r_moon_earth = -r_earth_moon;
    r_moon_sun = -r_sun_earth - r_earth_moon;
    r_moon_jupiter = r_sun_jupiter - r_sun_earth - r_earth_moon;
    fprintf('r_moon_earth = %.5e\n', r_moon_earth);
    fprintf('r_moon_sun = %.5e\n', r_moon_sun);
    fprintf('r_moon_jupiter = %.5e\n', r_moon_jupiter);

    % calculate dominant term for spacecraft relative to moon
    dominant = -(Gm_moon*r_moon_sc)/((abs(r_moon_sc))^3);
    fprintf('dominant = %.5e km/s^2\n\n', dominant);

    % due to earth
    direct_earth = Gm_earth*r_sc_earth/(abs(r_sc_earth))^3;
    indirect_earth = Gm_earth*r_moon_earth/(abs(r_moon_earth))^3;
    net_pert_earth = direct_earth - indirect_earth;
    fprintf('direct_earth = %.5e km/s^2\n', direct_earth);
    fprintf('indirect_earth = %.5e km/s^2\n', indirect_earth);
    fprintf('net_pert_earth = %.5e km/s^2\n\n', net_pert_earth);

    % due to sun
    direct_sun = Gm_sun*r_sc_sun/(abs(r_sc_sun))^3;
    indirect_sun = Gm_sun*r_moon_sun/(abs(r_moon_sun))^3;
    net_pert_sun = direct_sun - indirect_sun;
    fprintf('direct_sun = %.5e km/s^2\n', direct_sun);
    fprintf('indirect_sun = %.5e km/s^2\n', indirect_sun);
    fprintf('net_pert_sun = %.5e km/s^2\n\n', net_pert_sun);

    % due to jupiter
    direct_jupiter = Gm_jupiter*r_sc_jupiter/(abs(r_sc_jupiter))^3;
    indirect_jupiter = Gm_jupiter*r_moon_jupiter/(abs(r_moon_jupiter))^3;
    net_pert_jupiter = direct_jupiter - indirect_jupiter;
    fprintf('direct_jupiter = %.5e km/s^2\n', direct_jupiter);
    fprintf('indirect_jupiter = %.5e km/s^2\n', indirect_jupiter);
    fprintf('net_pert_jupiter = %.5e km/s^2\n\n', net_pert_jupiter);

    % total net acceleration
    ddot_r_moon_sc = dominant + net_pert_earth + net_pert_sun + net_pert_jupiter;
    fprintf('total net acceleration = %.5e km/s^2\n', ddot_r_moon_sc);
end


% spacecraft relative to earth
function earth_sc(Gm_earth, Gm_jupiter, Gm_moon, Gm_sun, ...
    r_moon_sc, r_sun_earth, r_earth_moon, r_sun_jupiter)

    % position vectors relative to spacecraft (sc)
    r_sc_moon = -r_moon_sc;
    r_sc_sun = -r_sun_earth - r_earth_moon - r_moon_sc;
    r_sc_jupiter = r_sun_jupiter - r_sun_earth - r_earth_moon - r_moon_sc;
    fprintf('r_sc_moon = %.5e\n', r_sc_moon);
    fprintf('r_sc_sun = %.5e\n', r_sc_sun);
    fprintf('r_sc_jupiter = %.5e\n', r_sc_jupiter);

    % position vectors relative to earth
    r_earth_sc = r_earth_moon + r_moon_sc;
    r_earth_moon = r_earth_moon;
    r_earth_sun = -r_sun_earth;
    r_earth_jupiter = r_sun_jupiter - r_sun_earth;
    fprintf('r_earth_sc = %.5e\n', r_earth_sc);
    fprintf('r_earth_moon = %.5e\n', r_earth_moon);
    fprintf('r_earth_sun = %.5e\n', r_earth_sun);
    fprintf('r_earth_jupiter = %.5e\n', r_earth_jupiter);

    % calculate dominant term for spacecraft relative to moon
    dominant = -(Gm_earth*r_earth_sc)/((abs(r_earth_sc))^3);
    fprintf('dominant = %.5e km/s^2\n\n', dominant);

    % due to moon
    direct_moon = Gm_moon*r_sc_moon/(abs(r_sc_moon))^3;
    indirect_moon = Gm_moon*r_earth_moon/(abs(r_earth_moon))^3;
    net_pert_moon = direct_moon - indirect_moon;
    fprintf('direct_moon = %.5e km/s^2\n', direct_moon);
    fprintf('indirect_moon = %.5e km/s^2\n', indirect_moon);
    fprintf('net_pert_moon = %.5e km/s^2\n\n', net_pert_moon);

    % due to sun
    direct_sun = Gm_sun*r_sc_sun/(abs(r_sc_sun))^3;
    indirect_sun = Gm_sun*r_earth_sun/(abs(r_earth_sun))^3;
    net_pert_sun = direct_sun - indirect_sun;
    fprintf('direct_sun = %.5e km/s^2\n', direct_sun);
    fprintf('indirect_sun = %.5e km/s^2\n', indirect_sun);
    fprintf('net_pert_sun = %.5e km/s^2\n\n', net_pert_sun);

    % due to jupiter
    direct_jupiter = Gm_jupiter*r_sc_jupiter/(abs(r_sc_jupiter))^3;
    indirect_jupiter = Gm_jupiter*r_earth_jupiter/(abs(r_earth_jupiter))^3;
    net_pert_jupiter = direct_jupiter - indirect_jupiter;
    fprintf('direct_jupiter = %.5e km/s^2\n', direct_jupiter);
    fprintf('indirect_jupiter = %.5e km/s^2\n', indirect_jupiter);
    fprintf('net_pert_jupiter = %.5e km/s^2\n\n', net_pert_jupiter);

    % total net acceleration
    ddot_r_earth_sc = dominant + net_pert_moon + net_pert_sun + net_pert_jupiter;
    fprintf('total net acceleration = %.5e km/s^2\n', ddot_r_earth_sc);
end