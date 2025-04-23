        if LATTEN then begin
          pTplHgtpatten:=FAPOLY*abs(Nucfit[Inuc].Coefficient)*TPL[Ipks].Eff*TplFactC;
          for J:=1 TO 6 do begin
            if (VaryTAttn0[j]) then begin
              Kpar:=Kpar+1;
              if (NUCfit[Inuc].Attenuation = 0)then Pfa[Kpar]:=Pfa[Kpar]+
               pTplHgtpatten*Tpl[Ipks].pAtten[J];
              end;
            end;
          for J:=1 TO 6 do begin
            if (VaryTAttn1[j]) then begin
              Kpar:=Kpar+1;
              if (NUCfit[Inuc].Attenuation = 1)then Pfa[Kpar]:=Pfa[Kpar]+
               pTplHgtpatten*Tpl[Ipks].pAtten[J];
              end;
            end;
          for J:=1 TO 6 do begin
            if (VaryTAttn2[j]) then begin
              Kpar:=Kpar+1;
              if (NUCfit[Inuc].Attenuation = 2)then Pfa[Kpar]:=Pfa[Kpar]+
               pTplHgtpatten*Tpl[Ipks].pAtten[J];
              end;
            end;
          end;

