VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)
LIGHT_SPEED = 299792458; % m/s

%% Plot water permitivity
% All functions using constants on p257 of Ellison et al. 257
% Note tau is in ps, fml
% Using middle set of parameters which indicate best fit

% clear; close all; clc

freqs0_100 = (0:0.5:100)*10^9;

temps = [10, 20, 25, 30, 40];
cols = jet(5);

% Param x temp x freq range
tempParam = cell(5, 5, 4);
% params are e0, e1, tau, lower freq, upper freq

% 10 deg
tempParam(:,1,1) = {83.928, 6.821, 12.81, 1, 9.5};
tempParam(:,1,2) = {83.134, 5.310, 12.41, 8.5, 20};
tempParam(:,1,3) = {88.316, 5.801, 13.38, 19, 100};
tempParam(:,1,4) = {NaN, NaN, NaN, NaN, NaN};

% 20
tempParam(:,2,1) = {80.181, 6.507, 9.63, 1, 9.5};
tempParam(:,2,2) = {79.989, 5.953, 9.54, 8.5, 20};
tempParam(:,2,3) = {80.467, 5.724, 9.45, 19, 40};
tempParam(:,2,4) = {83.055, 6.024, 9.69, 30, 100};

% 25
tempParam(:,3,1) = {78.404, 6.640, 8.45, 1, 9.5};
tempParam(:,3,2) = {77.657, 4.633, 8.14, 8.5, 20};
tempParam(:,3,3) = {78.746, 5.676, 8.38, 19, 40};
tempParam(:,3,4) = {78.814, 5.492, 8.26, 30, 100};

% 30
tempParam(:,4,1) = {76.582, 6.903, 7.56, 1, 9.5};
tempParam(:,4,2) = {76.452, 5.685, 7.38, 8.5, 20};
tempParam(:,4,3) = {76.867, 5.321, 7.18, 19, 100};
tempParam(:,4,4) = {NaN, NaN, NaN, NaN, NaN};

% 40
tempParam(:,5,1) = {73.290, 9.247, 6.32, 1, 9.5};
tempParam(:,5,2) = {72.919, 5.671, 5.96, 8.5, 100};
tempParam(:,5,3) = {NaN, NaN, NaN, NaN, NaN};
tempParam(:,5,4) = {NaN, NaN, NaN, NaN, NaN};

%% Plot real and imag permitivity
figure

for iTemp = 1:length(temps)
   for jRange = 1:4
       if ~isnan(tempParam{1, iTemp, jRange})
           e0 = tempParam{1, iTemp, jRange};      
           e1 = tempParam{2, iTemp, jRange};   
           tau = tempParam{3, iTemp, jRange}*10^-12;
           freqs = (tempParam{4, iTemp, jRange}:0.5:tempParam{5, iTemp, jRange})*10^9;
           
           permitivity = (e0-e1)./(1+1j*2*pi*freqs*tau)+e1;

           subplot(1,2,1); hold on; plot(freqs/10^9, real(permitivity), 'color', cols(iTemp,:))
           subplot(1,2,2); hold on; plot(freqs/10^9, -imag(permitivity), 'color', cols(iTemp,:))
       end
   end
end

subplot(1,2,1); title('Real')
subplot(1,2,2); title('Imaginary')

%% plot wavelength change
%%% Wavelength always long relative to sizes (drops/viruses)
figure; hold on
LIGHT_SPEED = 299792458; % m/s

glassPermitivity = 4; % not sure how constant this is across frequncy range

plot(freqs0_100/10^9, LIGHT_SPEED./freqs0_100, 'k')

plot(freqs0_100/10^9, LIGHT_SPEED./freqs0_100/sqrt(glassPermitivity), 'm')

for iTemp = 1:length(temps)
   for jRange = 1:4
       if ~isnan(tempParam{1, iTemp, jRange})
           e0 = tempParam{1, iTemp, jRange};      
           e1 = tempParam{2, iTemp, jRange};   
           tau = tempParam{3, iTemp, jRange}*10^-12;
           freqs = (tempParam{4, iTemp, jRange}:0.5:tempParam{5, iTemp, jRange})*10^9;
           
           permitivity = (e0-e1)./(1+1j*2*pi*freqs*tau)+e1;

           plot(freqs/10^9, LIGHT_SPEED./freqs./sqrt(real(permitivity)), 'color', cols(iTemp,:))
       end
   end
end

line([0 100], [1 1]*10^-3)
line([0 100], 100*[1 1]*10^-6)
line([0 100], [1 1]*10^-6)
line([0 100], 100*[1 1]*10^-9)
ylim([0 5*10^-3])

title('Wavelength in media')

%% Plot e-field in media

figure;

warning('calc of n in geometric should use complex modulus')

for iTemp = 1:length(temps)
   for jRange = 1:4
       if ~isnan(tempParam{1, iTemp, jRange})
           e0 = tempParam{1, iTemp, jRange};      
           e1 = tempParam{2, iTemp, jRange};   
           tau = tempParam{3, iTemp, jRange}*10^-12;
           freqs = (tempParam{4, iTemp, jRange}:0.5:tempParam{5, iTemp, jRange})*10^9;
           
           permitivity = real((e0-e1)./(1+1j*2*pi*freqs*tau)+e1);

           % in drop
           subplot(2,2,1); hold on; plot(freqs/10^9, 3./(permitivity+2), 'color', cols(iTemp,:))
           subplot(2,2,2); hold on; plot(freqs/10^9, 3*glassPermitivity./(permitivity+2*glassPermitivity), 'color', cols(iTemp,:))
       
           % geometric
           subplot(2,2,3); hold on; plot(freqs/10^9, (1-((1-sqrt(permitivity))./(1+sqrt(permitivity))).^2), 'color', cols(iTemp,:))
           subplot(2,2,4); hold on; plot(freqs/10^9, (1-((sqrt(glassPermitivity)-sqrt(permitivity))./ ...
               (sqrt(glassPermitivity)+sqrt(permitivity))).^2), 'color', cols(iTemp,:))
           
       end
   end
end

subplot(2,2,1); title('Water drop in air'); ylim([0 1])
subplot(2,2,2); title('Water drop in glass'); ylim([0 1])

subplot(2,2,3); title('Refraction from air into water'); ylim([0 1])
subplot(2,2,4); title('Refraction from glass into water'); ylim([0 1])

%% Plot power limits

figure; hold on;

plot(6:0.5:96, (200*((6:0.5:96)/3).^0.2)/2, 'b')
line([96 100], [1 1]*400/2, 'color', 'b')

line(6:0.5:30, (18.56*(6:0.5:30).^0.699)/2, 'color', 'g')
line([30 100], [1 1]*200/2, 'color', 'g')

line([6 100], [1 1]*50/2, 'color', 'b', 'linestyle', '--')
line([6 100], [1 1]*10/2, 'color', 'g', 'linestyle', '--')

legend('Safe power - peak, occupational', 'Safe power - peak, public', ...
    'Safe power - average, occupational', 'Safe power - average, public')

%% plot resonances
wetSpeed = 1920; % Longitudional speeds
drySpeed = 3430;

diameters = [80 100 120]*10^-9;

zerosToCalc = 12; % counting from zero

resonantFreqs = zeros(length(diameters) , zerosToCalc, 2);

for iDiam = 1:length(diameters)
    resonantFreqs(iDiam, :, 1) = calcualtesphereresonance(diameters(iDiam)/2, ...
                    'sph', 1, zerosToCalc-1, wetSpeed, wetSpeed/2, 5*10^9, 10^6, 0);
                
    resonantFreqs(iDiam, :, 2) = calcualtesphereresonance(diameters(iDiam)/2, ...
                    'sph', 1, zerosToCalc-1, drySpeed, drySpeed/2, 5*10^9, 10^6, 0);            
end

figure; hold on
for iZero = 1:zerosToCalc
    % plot wet
    line(resonantFreqs([1 3], iZero, 1)/10^9, [1 1]*0.5+iZero/zerosToCalc, 'color', 'm') 
    plot(resonantFreqs(2, iZero, 1)/10^9, 0.5+iZero/zerosToCalc, 'mo') 
    
    %plot dry
    line(resonantFreqs([1 3], iZero, 2)/10^9, [1 1]*1.5+iZero/zerosToCalc, 'color', 'k') 
    plot(resonantFreqs(2, iZero, 2)/10^9, 1.5+iZero/zerosToCalc, 'ko') 
end

xlim([0 100])
ylim([0 3])

%% Compare amplitude given shifting resonant frequencies
% Values from 13nm Sun's nanosphere, but frequncies change to ones above for nanosphere
reducedMass = 4.4494e-22;
systemQ = 2.0625e+00;
qToUse = 5.8171e-17;
nanocrystalNumber = 1.6000e+16;
apertureArea = 2.1000e-04;
   
amplitudeValues = zeros(zerosToCalc,1);
absorbtionValues = zeros(zerosToCalc,1);
powerValues = zeros(zerosToCalc,1);
extinctionValues = zeros(zerosToCalc,1);

for iZero = 1:zerosToCalc
    [absorbtionValues(iZero), extinctionValues(iZero), powerValues(iZero), amplitudeValues(iZero)] = calculatesphereabsorbtion(...
        resonantFreqs(2,iZero,1)*2*pi, resonantFreqs(2,iZero,1)*2*pi, ...
        reducedMass, systemQ, qToUse, nanocrystalNumber/2, apertureArea, 1, []);
end

figure; 
subplot(2,2,1)
plot(resonantFreqs(2,:,1)/10^9, amplitudeValues/max(amplitudeValues), 'x')
xlim([0 100]); ylim([0 1]); title('amplitude')

subplot(2,2,2)
plot(resonantFreqs(2,:,1)/10^9, powerValues/max(powerValues), 'x')
xlim([0 100]); ylim([0 1]); title('power')

subplot(2,2,3)
plot(resonantFreqs(2,:,1)/10^9, extinctionValues/max(extinctionValues), 'x')
xlim([0 100]); ylim([0 1]); title('extinction')

subplot(2,2,4)
plot(resonantFreqs(2,:,1)/10^9, absorbtionValues/max(absorbtionValues), 'x')
xlim([0 100]); ylim([0 1]); title('absoption')

%% Play with heat transfer

% Sun use 486 W/m^2 at 6 GHz for 15 minutes
% use model for 30 as went from 27 to 35 degress
e0 = tempParam{1, 4, 1};      
e1 = tempParam{2, 4, 1};   
tau = tempParam{3, 4, 1}*10^-12;
freq = 6*10^9;

permitivity = (e0-e1)./(1+1j*2*pi*freq*tau)+e1;

seconds = 15*60;
power = 486; % Should be average or peak?

% Equation 7 from Canals 1999 - for a droplet but size independent.
    %Peak rather than average amplitude seems to be used...
dTpS = 2*pi*freq*VACCUM_PERMITIVITY*imag(-permitivity)/997/4200*...
    sqrt(power/(LIGHT_SPEED*VACCUM_PERMITIVITY*sqrt(real(permitivity))))^2
%%% Why is square of field used here but power else where?

% Factor of x10 larger than Sun recorded - obviously there was some heat loss...
predictedHeat = dTpS*seconds
heatError = dTpS*seconds/7

power = 10;

%%% This has to be wrong - gives transmitted power larger than incident...
%%% In paper it is basically the same
2*pi*2.45*10^9*VACCUM_PERMITIVITY*9*800^2/((20.3*1000)^2*LIGHT_SPEED*VACCUM_PERMITIVITY)

% Try across freq and temp range
figure; hold on

for iTemp = 1:length(temps)
   for jRange = 1:4
       if ~isnan(tempParam{1, iTemp, jRange})
           e0 = tempParam{1, iTemp, jRange};      
           e1 = tempParam{2, iTemp, jRange};   
           tau = tempParam{3, iTemp, jRange}*10^-12;
           freqs = (tempParam{4, iTemp, jRange}:0.5:tempParam{5, iTemp, jRange})*10^9;
           
           permitivity = (e0-e1)./(1+1j*2*pi*freqs*tau)+e1;

           % eField for inside drop.
           eField = sqrt(power/(LIGHT_SPEED*VACCUM_PERMITIVITY)).*(3./(real(permitivity) + 2));
           
           dTpS = 2*pi.*freqs*VACCUM_PERMITIVITY.*imag(-permitivity)/997/4200.*eField.^2;

           subplot(1,2,1); hold on; plot(freqs/10^9, dTpS, 'color', cols(iTemp,:))
           subplot(1,2,2); hold on; plot(freqs/10^9, 1./(dTpS), 'color', cols(iTemp,:))
       end
   end
end

subplot(1,2,1); title('Predicted dt/ds (w/o error)')
subplot(1,2,2); title('Predicted ds for 1 dt (w/o error)')
ylim([0 10*60])

%% Play with mie absorbtion calculation
% Do for Douglas values first
dropletDiameter = [0.0005 0.1 0.5]*10^-2; % From Douglas 2004

%%% For bulk, try treating as 1mm radius to get absorbtion coeff.

figure;
for xDrop = 1:length(dropletDiameter)
    for iTemp = 1:length(temps)
       for jRange = 1
           if ~isnan(tempParam{1, iTemp, jRange})
               e0 = tempParam{1, iTemp, jRange};      
               e1 = tempParam{2, iTemp, jRange};   
               tau = tempParam{3, iTemp, jRange}*10^-12;
%                freqs = (tempParam{4, iTemp, jRange}:5:tempParam{5, iTemp, jRange})*10^9;
               freqs = [2.45 6 10]*10^9;
               
               for kFreq = 1:length(freqs)
                   
                   permitivity = (e0-e1)./(1+1j*2*pi*freqs(kFreq)*tau)+e1;
                    
                   % From: https://en.wikipedia.org/wiki/Refractive_index#Relative_permittivity_and_permeability
                   permitivity_mod = sqrt(real(permitivity)^2+imag(permitivity)^2);
                   
                   m = sqrt((permitivity_mod+real(permitivity))/2) +...
                       1j*sqrt((permitivity_mod-real(permitivity))/2);

                   x = pi*dropletDiameter(xDrop)./(LIGHT_SPEED/(freqs(kFreq)));

                   Qabs = mie_abs(m, x);
                   
                   extinction = Qabs*pi*dropletDiameter(xDrop)^2;
                   
                   subplot(6,2,xDrop*2-1); hold on; plot(freqs(kFreq)/10^9, Qabs, 'x', 'color', cols(iTemp,:))
                   
                   dTpS = extinction*power/997/4200;
                   
%                    subplot(6,2,xDrop*2); hold on; plot(freqs(kFreq)/10^9, dTpS, 'x', 'color', cols(iTemp,:))
                   %ylim([0 100*60])
               end
           end
       end
    end
end

% Max increase is around factor of 100 over 100 Ghz -> less than 10^6 diff in power...
    % Could have some effect, but probably not huge

% Do for all values now
% Using 10 mm here for 1ml drop size gives much smaller temp change than
% other eqn - 

dropletDiameter = [1 10 100 10000]*10^-6;

for xDrop = 1:length(dropletDiameter)
    for iTemp = 1:length(temps)
       for jRange = 1:4
           if ~isnan(tempParam{1, iTemp, jRange})
               e0 = tempParam{1, iTemp, jRange};      
               e1 = tempParam{2, iTemp, jRange};   
               tau = tempParam{3, iTemp, jRange}*10^-12;
               freqs = (tempParam{4, iTemp, jRange}:5:tempParam{5, iTemp, jRange})*10^9;
               
               for kFreq = 1:length(freqs)
                   
                   permitivity = (e0-e1)./(1+1j*2*pi*freqs(kFreq)*tau)+e1;
                    
                   % From: https://en.wikipedia.org/wiki/Refractive_index#Relative_permittivity_and_permeability
                   permitivity_mod = sqrt(real(permitivity)^2+imag(permitivity)^2);
                   
                   m = sqrt((permitivity_mod+real(permitivity))/2) +...
                       1j*sqrt((permitivity_mod-real(permitivity))/2);

                   x = pi*dropletDiameter(xDrop)./(LIGHT_SPEED/(freqs(kFreq)));

                   Qabs = mie_abs(m, x);
                   
                   extinction = Qabs*pi*dropletDiameter(xDrop)^2;
                   
                   subplot(6,2,xDrop*2); hold on; plot(freqs(kFreq)/10^9, Qabs, 'x', 'color', cols(iTemp,:))
                   
                   dTpS = extinction*power/997/4200;
                   
%                    subplot(6,2,xDrop*2); hold on; plot(freqs(kFreq)/10^9, dTpS, 'x', 'color', cols(iTemp,:))
                   %ylim([0 100*60])
               end
           end
       end
    end
end

%% Test for sphere embedded in material

dropletDiameter = (100:200:1000)*10^-6; %(10:10:200)*10^-6; % (100:200:1000)*10^-6;

colsD = jet(20);

glassPerm = 4.2 + 0.04j; % borosilicate glass
glass_mod = sqrt(real(glassPerm)^2+imag(glassPerm)^2);
                   
glassN = sqrt((glass_mod+real(glassPerm))/2) +...
   1j*sqrt((glass_mod-real(glassPerm))/2);
                   
acrylicPerm = 2.65 + 0.02j;
acrylic_mod = sqrt(real(acrylicPerm)^2+imag(acrylicPerm)^2);
                   
acrylicN = sqrt((acrylic_mod+real(acrylicPerm))/2) +...
   1j*sqrt((acrylic_mod-real(acrylicPerm))/2);

figure; 
subplot(2,3,1); hold on
plot(1:20, (LIGHT_SPEED./((1:20)*10^9)/real(glassN)))
plot(1:20, (LIGHT_SPEED./((1:20)*10^9)/real(acrylicN)))
plot(1:20, (LIGHT_SPEED./((1:20)*10^9)))

for xDrop = 1:length(dropletDiameter)
    subplot(2,3,4); hold on
    plot(1:20, (LIGHT_SPEED./((1:20)*10^9)/real(glassN))./dropletDiameter(xDrop))
    
    subplot(2,3,5); hold on
    plot(1:20, (LIGHT_SPEED./((1:20)*10^9)/real(acrylicN))./dropletDiameter(xDrop))
    
    subplot(2,3,6); hold on
    plot(1:20, (LIGHT_SPEED./((1:20)*10^9))./dropletDiameter(xDrop))    
end

for i = 1:3
   subplot(2,3,3+i); hold on
   ylim([0 100])
   line([1 20], [10 10])
   
end

figure

for xDrop = 1:length(dropletDiameter)
    for iTemp = 3 % just at 25 deg
       for jRange = 1:2 %Just to 20 GHz
           if ~isnan(tempParam{1, iTemp, jRange})
               e0 = tempParam{1, iTemp, jRange};      
               e1 = tempParam{2, iTemp, jRange};   
               tau = tempParam{3, iTemp, jRange}*10^-12;
               freqs = (tempParam{4, iTemp, jRange}:2:tempParam{5, iTemp, jRange})*10^9;
               
               for kFreq = 1:length(freqs)
                   
                   permitivity = (e0-e1)./(1+1j*2*pi*freqs(kFreq)*tau)+e1;
                    
                   permitivity_mod = sqrt(real(permitivity)^2+imag(permitivity)^2);
                   
                   waterN = sqrt((permitivity_mod+real(permitivity))/2) +...
                       1j*sqrt((permitivity_mod-real(permitivity))/2);

                   % For glass
                   m = waterN/glassN;
                   
                   x = pi*dropletDiameter(xDrop)./((LIGHT_SPEED/(freqs(kFreq)))/real(glassN));
                   
                   Qabs = mie_abs(m, x);
                   
                   subplot(3,3,1); hold on; plot(dropletDiameter(xDrop)*10^3, Qabs, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   nj=5*round(2+x+4*x.^(1/3))+160;
                   eValues = sqrt(mie_esquare(m, x, nj));
                   
                   % Take mean weighted by area
                   dx=x/nj; xj=(0:dx:x);
                   xarea = 4*pi*xj.^2;
                   
                   subplot(3,3,4); hold on; plot(dropletDiameter(xDrop)*10^3, wmean(eValues, xarea), 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   % limit should be correct if x << 1
                   subplot(3,3,7); hold on; plot(dropletDiameter(xDrop)*10^3, x, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))

                   
                   % For acrylic
                   m = waterN/acrylicN;
                   
                   x = pi*dropletDiameter(xDrop)./((LIGHT_SPEED/(freqs(kFreq)))/real(acrylicN));

                   Qabs = mie_abs(m, x);
                   
                   subplot(3,3,2); hold on; plot(dropletDiameter(xDrop)*10^3, Qabs, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   nj=5*round(2+x+4*x.^(1/3))+160;
                   eValues = sqrt(mie_esquare(m, x, nj));
                   
                   dx=x/nj; xj=(0:dx:x);
                   xarea = 4*pi*xj.^2;
                   
                   subplot(3,3,5); hold on; plot(dropletDiameter(xDrop)*10^3, wmean(eValues, xarea), 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   subplot(3,3,8); hold on; plot(dropletDiameter(xDrop)*10^3, x, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
  
                   
                   % For air
                   m = waterN;
                   
                   x = pi*dropletDiameter(xDrop)./((LIGHT_SPEED/(freqs(kFreq))));

                   Qabs = mie_abs(m, x);
                   
                   subplot(3,3,3); hold on; plot(dropletDiameter(xDrop)*10^3, Qabs, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   nj=5*round(2+x+4*x.^(1/3))+160;
                   eValues = sqrt(mie_esquare(m, x, nj));
                   
                   dx=x/nj; xj=(0:dx:x);
                   xarea = 4*pi*xj.^2;
                   
                   subplot(3,3,6); hold on; plot(dropletDiameter(xDrop)*10^3, wmean(eValues, xarea), 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
                   subplot(3,3,9); hold on; plot(dropletDiameter(xDrop)*10^3, x, 'x', 'color', colsD(round(freqs(kFreq)/10^9),:))
                   
               end
           end
       end
    end
end

for i = 1:3
   subplot(3,3,i);
   ylim([0 0.3])
   line([0 1], [0.1 0.1])
end

% These actualy vary as a function of freq because of change in perm.
% Using 10 Ghz and 25 o for now...
subplot(3,3,4);
line([0 1], [1 1]*(real(glassPerm)*3)/(real(glassPerm)*2+62.52))

subplot(3,3,5);
line([0 1], [1 1]*(real(acrylicPerm)*3)/(real(acrylicPerm)*2+62.52))

subplot(3,3,6);
line([0 1], [1 1]*(3)/(2+62.52))

%% Exp with heating of virus directly - bit of a mess now..

power = 10; %486

% 10% Attenuation at 6Ghz - would have come from subset of viruses...
ExAbs = -1/(1.25*10^-3*7.5*10^14)*log(1-0.1);

% From 1972 Structure book, hydrated
mass = 1.2*10^-18;

diameter = 120*10^-9;

% Reduced concentration to reflect lower number of smaller spheres
% dist of mass would also change, most have smaller mass than peak at 6GHz 
conc = 7.5*10^14*0.01;

deltaTm = 7;

%From image, looks like around 1ml
    % May be heating bulk, but almost certainly wouldn't heat drop - too few viruses
    % Unless density of viruses is very high.
vol = 1*10^-3; (4/3*pi*(10/2*10^-6)^3); 

nViruses = vol*conc;

kw = 0.6;

% Eqn 15 from Sci Rep paper
deltaTv = ExAbs*power/4/pi/(diameter/2)/kw

% From eqn 14
deltaTs = deltaTm/(nViruses);

Cwater = 4200;
Cvirus = 1300; % Big estimate

deltaTv2 = vol*Cwater/(mass*nViruses*Cvirus)*deltaTs

massSln = vol
virusMass = mass*nViruses

%%% Big diff between two deltaTvs suggests this wasn't at equilibirum

% Approx deltaTm is also crazy - would be unstable loss before that..
deltaTm2 = deltaTv/(massSln*Cwater)*virusMass*Cvirus*(nViruses)


%% Note sure what I was doing below here - somehting about heating

permitivity = 71.35 + 18.37j; % 30 deg at 6 GHz

airField = sqrt(power/(LIGHT_SPEED*VACCUM_PERMITIVITY));

dropField = 3*(real(permitivity)+2)*airField;

dropPower = dropField^2*LIGHT_SPEED*VACCUM_PERMITIVITY;

virusPowerAbs = ExAbs*dropPower

% From 50 micron water drop at 6 GHz
dropPowerAbs = 1.2*10^-5*(pi*(25*10^-6)^2)*power

% Not this is peak not avg field. no *2
field = sqrt(power/(LIGHT_SPEED*VACCUM_PERMITIVITY*sqrt(real(permitivity))));

PAbsVirusInd = sqrt(real(permitivity))*VACCUM_PERMITIVITY*LIGHT_SPEED*field^2*ExAbs;

PAbsVirusInd*7.5*10^14*(5*10^-6);

% Cross section of indvidual virus is x10 alrger than 50 micron water drop at 10 GHz
% 1.2*10^-5*(pi*(25*10^-6)^2)

%%% From gold NP paper
% sqrt(real(permitivity))*VACCUM_PERMITIVITY*LIGHT_SPEED*field^2*ExAbs/(4*pi*(50*10^-9)*0.6) 
% 5.1*10^4*10*10^-18/(4*pi*(15*10^-9)*0.6)

% just guessed density and specific heat
dtPs = PAbs/1100/3500

dT = dtPs*15*60

% Alt - they absorb 10% at density used in both abs and inactivation
PAbsVirusSol = sqrt(real(permitivity))*VACCUM_PERMITIVITY*LIGHT_SPEED*field^2*0.1

PAbsWater = pi*6*10^9*VACCUM_PERMITIVITY*imag(permitivity)*field^2

dtPsWater = PAbsWater/997/4200

dtWater = dtPsWater*15*60

% dtPs = PAbs/1100/3500
% 
% dT = dtPs*15*60