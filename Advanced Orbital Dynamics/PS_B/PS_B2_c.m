% For the characteristic length, we need the distance between the main body and the secondary body
% or semi-major axis of orbit will do
% for characteristic mass, we need the gravitational parameter of both, 
% divided by the gravitational constant
% for characteristic time, we need the char. length, and grav param of both
G = 6.67430e-20; % Constant of Gravitation kg^-1 km^3 s^-2

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
a_Charon = 17536; % about Pluto

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Jupiter = 126712767.8578;
mu_Mars = 42828.314258067;
mu_Saturn = 37940626.061137;
mu_Neptune = 6836534.0638793;
mu_Uranus = 5794549.0070719;
mu_Pluto = 981.600887707;
mu_Moon = 4902.8005821478;
mu_Phobos = 0.0007112;
mu_Titan = 8978.1382;
mu_Callisto = 7179.289;
mu_Oberon = 205.3;
mu_Triton = 1427.598;
mu_Charon = 102.30;

mass_table = table(zeros(0,2), cell(0,1), 'VariableNames', {'System', 'Char. Mass'});

m_Sun_Earth = (mu_Sun+mu_Earth)/G;
m_Sun_Earth_Data = {'m_Sun_Earth', m_Sun_Earth};
mass_table = [mass_table; m_Sun_Earth_Data]

m_Earth_Moon = (mu_Earth+mu_Moon)/G;
m_Earth_Moon_Data = {'m_Earth_Moon', m_Earth_Moon};
mass_table = [mass_table; m_Earth_Moon_Data]

m_Sun_Jupiter = (mu_Sun+mu_Jupiter)/G;
m_Sun_Jupiter_Data = {'m_Sun_Jupiter', m_Sun_Jupiter};
mass_table = [mass_table; m_Sun_Jupiter_Data]

m_Mars_Phobos = (mu_Mars+mu_Phobos)/G;
m_Mars_Phobos_Data = {'m_Mars_Phobos', m_Mars_Phobos};
mass_table = [mass_table; m_Mars_Phobos_Data]

m_Sun_Mars = (mu_Sun+mu_Mars)/G,
m_Sun_Mars_Data = {'m_Sun_Mars', m_Sun_Mars};
mass_table = [mass_table; m_Sun_Mars_Data]

m_Saturn_Titan = (mu_Saturn+mu_Titan)/G;
m_Saturn_Titan_Data = {'m_Saturn_Titan', m_Saturn_Titan};
mass_table = [mass_table; m_Saturn_Titan_Data]

m_Jupiter_Callisto = (mu_Jupiter+mu_Callisto)/G;
m_Jupiter_Callisto_Data = {'m_Jupiter_Callisto', m_Jupiter_Callisto};
mass_table = [mass_table; m_Jupiter_Callisto_Data]

m_Uranus_Oberon = (mu_Uranus+mu_Oberon)/G;
m_Uranus_Oberon_Data = {'m_Uranus_Oberon', m_Uranus_Oberon};
mass_table = [mass_table; m_Uranus_Oberon_Data]

m_Sun_Neptune = (mu_Sun+mu_Neptune)/G;
m_Sun_Neptune_Data = {'m_Sun_Neptune', m_Sun_Neptune};
mass_table = [mass_table; m_Sun_Neptune_Data]

m_Neptune_Triton = (mu_Neptune+mu_Triton)/G;
m_Neptune_Triton_Data = {'m_Neptune_Triton', m_Neptune_Triton};
mass_table = [mass_table; m_Neptune_Triton_Data]

m_Pluto_Charon = (mu_Pluto+mu_Charon)/G;
m_Pluto_Charon_Data = {'m_Pluto_Charon', m_Pluto_Charon};
mass_table = [mass_table; m_Pluto_Charon_Data]

disp(mass_table)

