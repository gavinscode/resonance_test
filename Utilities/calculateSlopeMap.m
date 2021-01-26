function [slopeMap, vL, vT] = calculateSlopeMap(vL, vT, modes, bigSize, smallSize)
 % Makes map of slopes across given range of Vl and Vt
 % can be used with single vL and vT, then just calculates slope for given
 % speeds
 
 % Velocity inputs should be in m/s, size inputs should be in m
 % Modes is +1 highest mode number as usual

    slopeMap = zeros(length(vL), length(vT), modes);

    for iLong = 1:length(vL)

        for jTrans = 1:length(vT)

            % 2nd and 3rd last parameters (guess and step size) are set to
            % work well for nanospheres, may change for viruses
            
            % Test two sizes to get slope (sizes aren't really important)
            highFreqs = calcualtesphereresonance(smallSize/2, ...
                'sph', 1, modes-1, vL(iLong), vT(jTrans), 5*10^9, 10^6, 0)*2*pi;

            lowFreqs = calcualtesphereresonance(bigSize/2, ...
                'sph', 1, modes-1, vL(iLong), vT(jTrans), 5*10^9, 10^6, 0)*2*pi;

            % Calc slope
            xTravel = 1/(smallSize/2) - 1/(bigSize/2);

            yTavel = highFreqs - lowFreqs;

            slopeMap(iLong, jTrans, :) = yTavel/xTravel;
        end
    end  

end

