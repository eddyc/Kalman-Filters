<CsoundSynthesizer>
<CsOptions>
-d
</CsOptions>
<CsInstruments>
sr = 44100
ksmps = 512
nchnls = 2
0dbfs = 1

ga_Phasor init 0

instr 1, Init

    a_Freq1 table ga_Phasor, 100, 1
    a_Freq2 table ga_Phasor, 101, 1
    a_Freq3 table ga_Phasor, 102, 1
    a_Freq4 table ga_Phasor, 103, 1
    
    k_Freq1 downsamp a_Freq1
    k_Freq2 downsamp a_Freq2
    k_Freq3 downsamp a_Freq3
    k_Freq4 downsamp a_Freq4

    a_Out1 vco2 0.2, k_Freq1
    a_Out2 vco2 0.2, k_Freq2
    a_Out3 vco2 0.2, k_Freq3
    a_Out4 vco2 0.2, k_Freq4

    a_Out = (a_Out1 + a_Out2 + a_Out3 + a_Out4) / 4
    
    a_Out gauss 1
    
    if  k_Freq1 < 0 then
    
        ;k_Freq1 = 0
        
    endif

    a_Res1 butterbp a_Out, k_Freq1, 20
    a_Res2 butterbp a_Out, k_Freq2, 20
    a_Res3 butterbp a_Out, k_Freq3, 20
    a_Res4 butterbp a_Out, k_Freq4, 20

    a_Out = (a_Res1 + a_Res2 + a_Res3 + a_Res4) 

    outs a_Out, a_Out
endin

instr 2, Phasor

    i_Increment = p4
    ga_Phasor phasor sr / ksmps / i_Increment

endin


</CsInstruments>  
<CsScore>
f 1 0 16384 10 1
e3600
</CsScore>
</CsoundSynthesizer>
