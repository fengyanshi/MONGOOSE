        parameter (mbub=10,Rgas=8.314)
        common /bub1/ dcl, dcr, dct, dcb
     &               ,rncl, rncr, rnct, rncb
        common /bub2/ cbub(mbub,nxy),bcprod(mbub,nxy)
     &         , bcdedu(mbub,nxy),bnutt(mbub,nxy),buoy(nxy)
     &         , rnum_b(mbub,nxy),bnprod(mbub,nxy),bndedu(mbub,nxy)
        common /bub3/ wb(mbub),rb(mbub),alpha_g(nxy),Temp
     &         ,rb_update(mbub,nxy),cbub_transf(mbub,nxy)
     &         ,rnum_transf(mbub,nxy),prod_b(mbub)
     &         ,prod_b_per_xprod_dt(mbub)
       

