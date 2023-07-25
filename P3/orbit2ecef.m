function [coordinate_eci_navigation,coordinate_ecef_navigation  ]= orbit2ecef(data)
%______________________________________________________________________________________________________%
%%%  this function perduced by AmirAbbas Saberi %%%
%----- input orbit and kepler parameter 
%----- to convert data from ORBITAL TO ECEF 
%----- output ECEF Coordinate data and ECI coordinate data
%----- output CalenderTime ----------------------------------
%______________________________________________________________________________________________________%
%%==== Global Constant =====%%
for t1 = 1:numel(data(:,1))/8

        prn = data(-7+8*t1,1);
        crs = data(-6+8*t1,2); 
        delta_n = data(-6+8*t1,3);
        M0 = data(-6+8*t1,4);
        cuc =  data(-5+8*t1,1); 
        e = data(-5+8*t1,2); 
        cus = data(-5+8*t1,3); 
        sqrtA = data(-5+8*t1,4);
        toe = data(-4+8*t1,1); 
        cic = data(-4+8*t1,2); 
        OMEGA =data(-4+8*t1,3);
        cis = data(-4+8*t1,4);
        i0 = data(-3+8*t1,1); 

        omega = data(-3+8*t1,3);
        OMEGA_DOT = data(-3+8*t1,4);
        i_dot = data(-2+8*t1,1) ;
        crc = data(-3+8*t1,2);

        year = data(-7+8*t1,2) ;
        month = data(-7+8*t1,3) ;
        day = data(-7+8*t1,4) ;
        hour = data(-7+8*t1,5) ;
        min = data(-7+8*t1,6) ;
        sec = data(-7+8*t1,7) ;
        
        JD = juliand(year,month,day,hour,min,sec) ;
        JD0 = juliand(1980,1,6,0,0,0) ;
        gps_week = floor( (JD - JD0) / 7 );

        gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend

        t = gps_seconds - toe;
        if t<0
            t = t+24*60*60*7;
        end

        GM = 3.986004418e+14;
        we = 7.2921151467e-5;
        a = sqrtA^2;


        M = M0 + (sqrt(GM/a^3) + delta_n)*t;

        if or((-pi < M < 0) , (M > pi))
            E = M - e;
        else
            E = M + e;
        end
        E_new=1000;
        while abs(E_new - E) > 10^-15
            E_new = E + ((M - E + e .* sin(E))./(1 - e .* cos(E)));
  
            E = E_new;
        end

        v = atan2((sqrt(1-e^2)*sin(E)),(cos(E) - e));

        u = omega +v+ cuc*cos(2*(omega + v)) + cus*sin(2*(omega + v));
        r = a*(1-e*cos(E)) + crc*cos(2*(omega + v)) + crs*sin(2*(omega + v));
        i = i0 + i_dot*t + cic*cos(2*(omega + v)) + cis*sin(2*(omega + v));
        omg = OMEGA + (OMEGA_DOT - we)*t - we*toe;

        v_pf = [r; 0; 0];

        R_ecef =[(cos(omg)*cos(u))-(sin(omg)*sin(u)*cos(i)) , (-cos(omg)*sin(u))-(sin(omg)*cos(u)*cos(i)) , (sin(omg)*sin(i));
         (sin(omg)*cos(u))+(cos(omg)*sin(u)*cos(i)) , (-sin(omg)*sin(u))+(cos(omg)*cos(u)*cos(i)) , (-cos(omg)*sin(i));
          sin(u)*sin(i)  ,  cos(u)*sin(i)   ,  cos(i)];
        
        ecef = R_ecef*v_pf;
        x = ecef(1,1);           
        y = ecef(2,1);
        z = ecef(3,1);
        coordinate_ecef_navigation(t1,1:12) = [prn,x,y,z,year,month,day,hour,min,sec,gps_seconds,toe]; %% Navigation ECEF Coordinate 

        R_eci= [cos(OMEGA) sin(OMEGA) 0;-sin(OMEGA) cos(OMEGA) 0;0 0 1]...
                  *[1 0 0;0 cos(i0) sin(i0);0 -sin(i0) cos(i0)]...
                  *[cos(omega) sin(omega) 0;-sin(omega) cos(omega) 0;0 0 1];
        
        eci = R_eci*[r*cos(v) ;r*sin(v); 0];
        x = eci(1,1);           
        y = eci(2,1);
        z = eci(3,1); 
        coordinate_eci_navigation(t1,1:12) = [prn,x,y,z,year,month,day,hour,min,sec,gps_seconds,toe]; %% Navigation ECEF Coordinate 


end
