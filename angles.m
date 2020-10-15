function Oz = angles(flag, n)
%Produces the angles over which we perform evaluations. Two possibilities:
%flag = 'discrete' gives Oz = [-1 1] and flag = 'lsgq' gives the S_n level 
%set gauss quadrature points for n = 2, 4, 6, 8, 10, 12, 14, 16, 18, 20.
if strcmpi(flag, 'discrete') == 1
    Oz = [-1 1];
    
elseif strcmpi(flag, 'lsgq') == 1
    Oz = lsgq(n);
       
else
    error('Angle option not available')   
end


    function angle = lsgq(k) 
       S2 = [-1 1];
       S4 = [-.3500212  -.8688903 .3500212 .8688903];
       S6 = [-.2666354 -.6815077 -.9261809 .266354 .6815077 .9261809];
       
       S8_p = [.2182179 .5773503 .7867958 .951189] ;
       S8_m = - S8_p;
       S8 = [S8_m S8_p];

      S10_p = [.18932213 .5088818 .7867958 .8397600 .9634910] ;
      S10_m = - S10_p;
      S10 = [S10_m S10_p];
      
      S12p = [.1672127 .4595476 .6280191 .7600210 .8722705 .9716377];
      S12m = - S12p;
      S12 = [S12m S12p];
      
      S14p = [.1519859 .4221570 .5773503 .6988921 .8022263 .8936911...
          .9766272];
      S14m = -S14p;
      S14 = [S14m S14p];
      
      S16p = [.138959 .3922893 .5370966 .6504265 .7467506 .8319966...
           .9092855 .9805009];
      S16m = - S16p;
      S16 = [S16m S16p];
      
      S18p = [.1293445 .3680438 .5041652 .6106625 .7011669 .7812562...
          .8538662  .9207680 .8931277];
      S18m = -S18p;
      S18 = [S18m S18p];
      
      S20p = [.1206033 .3475743 .4765193 .5773503 .6630204 .7388226...
          .8075404 .8708526 .9298639 .9853475];
      S20m = -S20p;
      S20 = [S20m S20p]; 
      
      Sn = {S2, S4, S6, S8, S10, S12, S14, S16, S18, S20};
      angle = Sn{k/2};
      
    end




end