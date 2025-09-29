% For the characteristic length, we need the distance between the main body and the secondary body
% or semi-major axis of orbit will do
% for characteristic mass, we need the gravitational parameter of both, 
% divided by the gravitational constant
% for characteristic time, we need the char. length, and grav param of both
G = 6.67430e-20; % Constant of Gravitation kg^-1 km^3 s^-2
AU = 149597870.7; % km
% Characteristic length
a_Earth = 149597898;
a_Jupiter = 778279959;
a_Mars = 227944135;
a_Neptune = 4498337290;
a_Moon = 384400; % around Earth
a_Phobos = 9376; % about Mars
a_Titan = 1221865; % about Saturn
a_Callisto = 1882700; % about Jupiter
a_Oberon = 583511; % about Uranus [https://ssd.jpl.nasa.gov/sats/elem/]
a_Triton = 354759; % about Neptune
% a_Charon = 17536; % about Pluto
a_Charon = 19600;

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Jupiter = 126712767.8578;
mu_Mars = 42828.314258067;
mu_Saturn = 37940626.061137;
mu_Neptune = 6836534.0638793;
mu_Uranus = 5794549.0070719;
% mu_Pluto = 981.600887707;
mu_Pluto = 13024.6e18*G;
mu_Moon = 4902.8005821478;
mu_Phobos = 0.0007112;
mu_Titan = 8978.1382;
mu_Callisto = 7179.289;
mu_Oberon = 205.3;
mu_Triton = 1427.598;
mu_Charon = 106.1;

% Orbital Periods for Checking [days]
P_Earth = 31558205/86400;
P_Moon = 2360592/86400;
P_Jupiter = 374479305/86400;
P_Phobos = 27552.96/86400;
P_Mars = 59356281/86400;
P_Titan = 1377648/86400;
P_Callisto = 1441931.19/86400;
P_Oberon = 13.463237; % given in days
P_Neptune = 5203578080/86400;
P_Triton = 5.876994; % given in days
P_Charon = 551856.7066/86400;

mass_table = cell2table(cell(0,7), 'VariableNames', {'System', 'Chars. Length', 'Char. Mass', 'mu', 'char. time', 'Period (check)', 'Period (Data)'});


l_Sun_Earth = a_Earth;
m_Sun_Earth = (mu_Sun+mu_Earth)/G;
mu_Sun_Earth = mu_Earth/(mu_Sun+mu_Earth);
t_Sun_Earth = sqrt(a_Earth^3/(mu_Sun+mu_Earth));
m_Sun_Earth_Data = {'m_Sun_Earth', a_Earth, m_Sun_Earth, mu_Sun_Earth, t_Sun_Earth/86400, t_Sun_Earth/86400*2*pi, P_Earth};
mass_table = [mass_table; m_Sun_Earth_Data];

l_Earth_Moon = a_Moon;
m_Earth_Moon = (mu_Earth+mu_Moon)/G;
mu_Earth_Moon = mu_Moon/(mu_Earth+mu_Moon);
t_Earth_Moon = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
m_Earth_Moon_Data = {'Earth_Moon', a_Moon, m_Earth_Moon, mu_Earth_Moon, t_Earth_Moon/86400, t_Earth_Moon/86400*2*pi, P_Moon};
mass_table = [mass_table; m_Earth_Moon_Data];

l_Sun_Jupiter = a_Jupiter/AU;
m_Sun_Jupiter = (mu_Sun+mu_Jupiter)/G;
mu_Sun_Jupiter = mu_Jupiter/(mu_Sun+mu_Jupiter);
t_Sun_Jupiter = sqrt(a_Jupiter^3/(mu_Sun+mu_Jupiter));% in years
m_Sun_Jupiter_Data = {'Sun_Jupiter', l_Sun_Jupiter, m_Sun_Jupiter, mu_Sun_Jupiter, t_Sun_Jupiter/86400/365, t_Sun_Jupiter/86400*2*pi, P_Jupiter};
mass_table = [mass_table; m_Sun_Jupiter_Data];

l_Mars_Phobos = a_Phobos/a_Moon; % a_Moon
m_Mars_Phobos = (mu_Mars+mu_Phobos)/G;
mu_Mars_Phobos = mu_Phobos/(mu_Mars+mu_Phobos);
t_Mars_Phobos = sqrt(a_Phobos^3/(mu_Mars+mu_Phobos)); % in hours
m_Mars_Phobos_Data = {'Mars_Phobos', l_Mars_Phobos, m_Mars_Phobos, mu_Mars_Phobos, t_Mars_Phobos/3600, t_Mars_Phobos/86400*2*pi, P_Phobos};
mass_table = [mass_table; m_Mars_Phobos_Data];

l_Sun_Mars = a_Mars/AU; % AU
m_Sun_Mars = (mu_Sun+mu_Mars)/G;
mu_Sun_Mars = mu_Mars/(mu_Sun+mu_Mars);
t_Sun_Mars = sqrt(a_Mars^3/(mu_Sun+mu_Mars));
m_Sun_Mars_Data = {'Sun_Mars', l_Sun_Mars, m_Sun_Mars, mu_Sun_Mars, t_Sun_Mars/86400, t_Sun_Mars/86400*2*pi, P_Mars};
mass_table = [mass_table; m_Sun_Mars_Data];

l_Saturn_Titan = a_Titan/a_Moon; % a_moon
m_Saturn_Titan = (mu_Saturn+mu_Titan)/G;
mu_Saturn_Titan = mu_Titan/(mu_Saturn+mu_Titan);
t_Saturn_Titan = sqrt(a_Titan^3/(mu_Saturn+mu_Titan));
m_Saturn_Titan_Data = {'Saturn_Titan', l_Saturn_Titan, m_Saturn_Titan, mu_Saturn_Titan, t_Saturn_Titan/86400, t_Saturn_Titan/86400*2*pi, P_Titan};
mass_table = [mass_table; m_Saturn_Titan_Data];

l_Jupiter_Callisto = a_Callisto/a_Moon;% a_moon
m_Jupiter_Callisto = (mu_Jupiter+mu_Callisto)/G;
mu_Jupiter_Callisto = mu_Callisto/(mu_Jupiter+mu_Callisto);
t_Jupiter_Callisto = sqrt(a_Callisto^3/(mu_Jupiter+mu_Callisto)); 
m_Jupiter_Callisto_Data = {'Jupiter_Callisto', l_Jupiter_Callisto, m_Jupiter_Callisto, mu_Jupiter_Callisto, t_Jupiter_Callisto/86400, t_Jupiter_Callisto/86400*2*pi, P_Callisto};
mass_table = [mass_table; m_Jupiter_Callisto_Data];

l_Uranus_Oberon = a_Oberon/a_Moon;% a_moon
m_Uranus_Oberon = (mu_Uranus+mu_Oberon)/G;
mu_Uranus_Oberon = mu_Oberon/(mu_Uranus+mu_Oberon);
t_Uranus_Oberon = sqrt(a_Oberon^3/(mu_Uranus+mu_Oberon));
m_Uranus_Oberon_Data = {'Uranus_Oberon', l_Uranus_Oberon, m_Uranus_Oberon, mu_Uranus_Oberon, t_Uranus_Oberon/86400, t_Uranus_Oberon/86400*2*pi, P_Oberon};
mass_table = [mass_table; m_Uranus_Oberon_Data];

l_Sun_Neptune = a_Neptune/AU;
m_Sun_Neptune = (mu_Sun+mu_Neptune)/G;
mu_Sun_Neptune = mu_Neptune/(mu_Sun+mu_Neptune);
t_Sun_Neptune = sqrt(a_Neptune^3/(mu_Sun+mu_Neptune)); % in years
m_Sun_Neptune_Data = {'Sun_Neptune', l_Sun_Neptune, m_Sun_Neptune, mu_Sun_Neptune, t_Sun_Neptune/86400/365, t_Sun_Neptune/86400*2*pi, P_Neptune};
mass_table = [mass_table; m_Sun_Neptune_Data];

l_Neptune_Triton = a_Triton/a_Moon; % a_moon
m_Neptune_Triton = (mu_Neptune+mu_Triton)/G;
mu_Neptune_Triton = mu_Triton/(mu_Neptune+mu_Triton);
t_Neptune_Triton = sqrt(a_Triton^3/(mu_Neptune+mu_Triton)); % in hours
m_Neptune_Triton_Data = {'Neptune_Triton', l_Neptune_Triton, m_Neptune_Triton, mu_Neptune_Triton, t_Neptune_Triton/3600, t_Neptune_Triton/86400*2*pi, P_Triton};
mass_table = [mass_table; m_Neptune_Triton_Data];

l_Pluto_Charon = a_Charon/a_Moon; % a_moon
m_Pluto_Charon = (mu_Pluto+mu_Charon)/G;
mu_Pluto_Charon = mu_Charon/(mu_Pluto+mu_Charon);
t_Pluto_Charon = sqrt(a_Charon^3/(mu_Pluto+mu_Charon));
m_Pluto_Charon_Data = {'Pluto_Charon', l_Pluto_Charon, m_Pluto_Charon, mu_Pluto_Charon, t_Pluto_Charon/86400, t_Pluto_Charon/86400*2*pi, P_Charon};
mass_table = [mass_table; m_Pluto_Charon_Data];

% format shortE
format short
disp(mass_table)

