function [bzj,czj] = func_Jpole(J_opt)
% @Description: J-pole expansion (from Xie2016, PST).
% @Filename: func_Jpole.m
% @Date: 2021-10-27
% @LastEditors: Bai Wei
% @LastEditTime: 2024-10-16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J=2.1:  Huba2009, 2-pole
% J=2.2:  Martin1979, 2-pole, J=2, I=3
% J=3:    Martin1980, 3-pole, J=3, I=3
% J=4.1:  Martin1980, 4-pole, J=4, I=5 
% J=4.2:  new calculation, J=4, I=5
% J=8.1:  Ronnmark1982, 8-pole for Z function
% J=8.2:  J=8, I=8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = floor(J_opt);

if(J==2)
    % opt=2; 
    if (J_opt==2.1)
        % Huba2009, 2-pole
        bzj(1) = -(0.5+0.81i);
        czj(1) = 0.51-0.81i;
    elseif (J_opt==2.2)
        % Martin1979, 2-pole, J=2, I=3
        bzj(1)=-(0.5+1.2891i);
        czj(1)=0.5138-1.0324i;
    else
        disp(['The opt=', num2str(round(mod(J_opt, 1)*10)),... 
              ' for J=', num2str(floor(J_opt)), ' does not exist, please reset!']);
    end
    bzj(2)=conj(bzj(1));
    czj(2)=-conj(czj(1));
    
elseif(J==3)
    % Martin1980, 3-pole, J=3, I=3
    bzj(1) = 0.1822+0.5756i;
    bzj(2) = -1.3643;
    czj(1) = -0.9217-0.9091i;
    czj(2) = -1.0204i;
    bzj(3) = conj(bzj(1));
    czj(3) = -conj(czj(1));

% J=4; % J-pole
elseif (J==4)
    % opt=2;
    if(J_opt==4.1)
        % Martin1980, 4-pole, J=4, I=5
        bzj(1) = 0.5468-0.0372i;
        bzj(2) = -1.0468+2.1018i;
        czj(1) = -1.2359-1.2150i;
        czj(2) = -0.3786-1.3509i;
    elseif(J_opt==4.2) % new calculation, J=4, I=5
        bzj(1) = 0.546796859834032 + 0.037196505239277i;
        bzj(2) = -1.046796859834027 + 2.101852568038518i;
        czj(1) = 1.23588765343592 - 1.21498213255731i;
        czj(2) = -0.378611612386277 - 1.350943585432730i;
    else
        disp(['The opt=', num2str(round(mod(J_opt, 1)*10)),... 
              ' for J=', num2str(floor(J_opt)), ' does not exist, please reset!']);
    end
    bzj(3:4) = conj(bzj(1:2));
    czj(3:4) = -conj(czj(1:2));
    
elseif (J==8)    
    % J-pole, usually J=8 is sufficient; other choice: 4, 12
    % opt=4;
    if(J_opt==8.1)
        % Ronnmark1982, 8-pole for Z function, and Z', J=8, I=10
        bzj(1) = -1.734012457471826E-2-4.630639291680322E-2i;
        bzj(2) = -7.399169923225014E-1+8.395179978099844E-1i;
        bzj(3) = 5.840628642184073+9.536009057643667E-1i;
        bzj(4) = -5.583371525286853-1.120854319126599E1i;
        czj(1) = 2.237687789201900-1.625940856173727i;
        czj(2) = 1.465234126106004-1.789620129162444i;
        czj(3) = 0.8392539817232638-1.891995045765206i;
        czj(4) = 0.2739362226285564-1.941786875844713i;
    elseif(J_opt==8.2) % J=8, I=8
        bzj(1)=  2.262756990456044 + 0.574863654552974i;
        bzj(2)=  0.007755112891899 + 0.435320280255620i;
        bzj(3)=  -0.017631167787210 - 0.001036701965463i;
        bzj(4)=  -2.752880935560734 - 4.080667522419399i;
        czj(1)=  -0.845779591068908 - 1.626106289905113i;
        czj(2)=  1.485551642694134 - 1.512332910786984i;
        czj(3)=  2.286671582218345 - 1.340185584073022i;
        czj(4)=  0.275217303800922 - 1.682797352333194i;
    else
        disp(['The opt=', num2str(round(mod(J_opt, 1)*10)),... 
              ' for J=', num2str(floor(J_opt)), ' does not exist, please reset!']);
    end
    bzj(5:8) = conj(bzj(1:4));
    czj(5:8) = -conj(czj(1:4));
else
    disp(['The ', 'J=', num2str(J), ' does not exist, please reset!']);
end

end
