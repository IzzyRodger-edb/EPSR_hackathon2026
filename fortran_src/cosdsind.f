      function cosd(deg)
      parameter (degradconv=3.141592653/180.0)
      cosd=cos(deg*degradconv)
      return
      end
      function sind(deg)
      parameter (degradconv=3.141592653/180.0)
      sind=sin(deg*degradconv)
      return
      end
      function acosd(cosphi)
      parameter (degradconv=3.141592653/180.0)
      acosd=acos(cosphi)
      acosd=acosd/degradconv
      return
      end   
      function asind(sinphi)
      parameter (degradconv=3.141592653/180.0)
      asind=asin(sinphi)
      asind=asind/degradconv
      return
      end   
      function atan2d(sinphi,cosphi)
      parameter (degradconv=3.141592653/180.0)
      atan2d=atan2(sinphi,cosphi)
      atan2d=atan2d/degradconv
      return
      end   