clear; 

nForceTol = 100; % Newtons, beginning force for 
% selection of end of first impact.
myrho = 0.05;
pass_frequency = 2500;
sample_rate = 44100;
dist_b4 = 0;
dist_af = 100;

%Data sets to look at
first = 1
last = 19
%Array for averaging the energy and velocity
eArray = zeros((last-first),1);
vArray = zeros((last-first),1);
i=1; %index to fill eArray and vArray

for cn = first:last %
%Read in data after running thru SensorData.py
    f = table2array(readtable("FF"+string(cn-1)+".csv"));
    a = table2array(readtable("AA"+string(cn-1)+".csv"));
    ID_table = table2array(readtable("II"+string(cn-1)+".csv"));

%Initialize Vars
    t= (0:(length(f)-2)).*1/44100;
    v = zeros([length(a),1]);

%Unit Convert to SI
    f = f.*4.44822; % lbf to Newtons
    a = a.*0.0254; % in/s^2 to m/s^2   

%Filter Data
    f_lowpass = lowpass(f,pass_frequency,sample_rate);
    f=f_lowpass;
    a_denoise = denoisetv(a,myrho);
    a_lowpass = lowpass(a,pass_frequency,sample_rate);
    a = a_lowpass;
    
%Rename some vars
    ID = ID_table(1);
    CID = ID_table(2);
    AID = ID_table(3);

    %Find the minimum acceleration index to find the moment of impact
    [~, minAccelIndex] = min(a(ID:CID));
    minAccelIndex=minAccelIndex+ID; %Add the ID value to offset correctly
    velocityAtImpact = trapz(a(ID:minAccelIndex)).*1/44100; %Integral from beginning to the moment of impact m/s
    vArray(i) = velocityAtImpact; %Store value for averaging
    fprintf("The velocity at impact is %f m/s\n ",velocityAtImpact)

    %Cumulative integral for the instant velo value
    cumVelo=cumtrapz(a(ID:CID));

    %Subplots for debugging
    %subplot(2,1,1)
    %plot(cumVelo)
    %title("Cumulative Velocity")
    %subplot(2,1,2)
    %plot(a(ID:CID+200))
    %title("Acceleration")
    %pause()

    %Calculate the cumulative velocity of impact
    v((ID):CID) = cumtrapz(a((ID):CID)).*1/44100;
    %fprintf('ID: %f\n',ID);
    %fprintf('CID: %f\n', CID);

    %Integrate for position
    s = cumtrapz(v).*1/44100./.0254;
    
    %Set any force less than 0 to 0?
    f(f<0) = 0;
    %Calculate force time velocity
    fv = f.*v;
    
    %Grabbing 120 data points?
    fvA = fv(ID:ID+120);
    fA = f(ID:ID+120);
    
    % Logical to isolate areas of tensile force -- SS
    % Areas of tensile force would not contribute any energy -- SS
    %Boolean vector? with 1 for if fv value is 0 and force is 0? -- JMC
    zr = (fvA == 0) & (fA == 0);
    % Calculate the derivative of FV to find minima maxima -- SS
    fvdiff = [0;diff(fv(ID:CID))];
    %Initialize some vars
    zcn = 0;
    nFTolFlag = 0;
    upflag = 0;
    for idc = 15:(length(zr)-1)
        % Set flag if the initial force tolerance is reached -- SS
       if fA(idc) >= nForceTol
           nFTolFlag = 1;
       end
       % These two ifs just find a minimum locaion -- SS
       if (fvA(idc) < -1) && ((fvA(idc)-fvA(idc+1))<0)
          upflag = 1; 
       end
       if upflag && ((fvA(idc)-fvA(idc+1))>=0)
          minloc = idc;
          break;
       end

       if zr(idc) && nFTolFlag && (fvdiff(idc) > 0 || fvdiff(idc-1) > 0 || fvdiff(idc+1) > 0)
           minloc = idc;
           break;
       end
    end
    %Subtracting "minloc" to only get plastic deformation energy?
    CID2 = ID+minloc;

    %Eliminate weird back samples
    for scn = 1:10
       if (fvA(scn) - fvA(scn+1)) > 0
          ID = ID + 1; 
       end
    end
    %Section to find min after the first small hump
    
    %Cumulative energy
    e = cumtrapz(fv(ID:CID)).*1/44100; %J
    %Total energy of impact (elastic and plastic)
    eOfImp=trapz(fv(ID:minAccelIndex)).*1/44100 %J
    fprintf("The energy of impact is %f J\n", eOfImp)
    eArray(i) = eOfImp;
    i=i+1;
    
%Plot everything
    subplot(4,1,1)
    hold on
    plot(fv(ID:CID))
    title('F*V Plot (Nm/s or Watts)')
    hold off
    subplot(4,1,2)
    plot(f(ID:CID))
    title('Force (N)')
    subplot(4,1,3)
    plot(a(ID:CID))
    title('Acceleration (m/s^2)')
    subplot(4,1,4)
    plot(v(ID:CID))
    title('Velocity (m/s)')
    %pause()

    %Velo plot
    %subplot(4,1,1)
    %plot(v(ID-50:CID))
    %title('Velocity(m/s)')
    %subplot(4,1,2)
    %plot(a(ID-50:CID))
    %plot(ID,1,"*")
    %pause()

end

%Average energy of impact
eArray = eArray * -1000
avgEOfImp = (sum(eArray)/length(eArray));
avgVOfImp = (sum(vArray)/length(vArray))*-1;
stdDevOfE = std(eArray);
stdDevOfV = std(vArray);
T = table(avgEOfImp,stdDevOfE,avgVOfImp,stdDevOfV);
T.Properties.VariableNames = {'Avg. Energy of Imp. (mJ)', 'Std. Dev. of Energy (mJ)', 'Average Velo. at Imp. (m/s)', 'Std. Dev. of Velo. (m/s)'}

%Conversions
%m/s to fps
%3.28

%Steel ball equilavent
densityOfBall = 0.00000793; %kg/mm^3
volOfBall = (4/3)*pi()*6.35^3; %mm^3
massofBall = densityOfBall*volOfBall; %kg
veloOfBall = sqrt((2*avgEOfImp)/massofBall); %m/s
fprintf("An equilavent energy of impact is a 304 SS ball traveling at %f m/s \n", veloOfBall)



%Function declarations
function f = denoisetv(g,mu)
I = length(g);
u = zeros(I,1);
y = zeros(I,1);
rho = 10;

eigD = abs(fftn([-1;1],[I 1])).^2;
for k=1:100
    f = real(ifft(fft(mu*g+rho*Dt(u)-Dt(y))./(mu+rho*eigD)));
    v = D(f)+(1/rho)*y;
    u = max(abs(v)-1/rho,0).*sign(v);
    y = y - rho*(u-D(f));
end
end

function y = D(x)
y = [diff(x);x(1)-x(end)];
end

function y = Dt(x)
y = [x(end)-x(1);-diff(x)];
end