using Base.Test
@test roundste(100.,34567.) == "100 ± 34600 (se)"
@test roundste(0.,12345.) == "0 ± 12300 (se)"
@test roundste(100.,4567.) == "100 ± 4570 (se)"
@test roundste(1099.,4567.) == "1100 ± 4570 (se)"
@test roundste(1049.,45649.) == "1000 ± 45600 (se)"
@test roundste(1049.,0.001) == "1049.00000 ± 0.00100 (se)"
@test roundste(1049.,0.1) == "1049.000 ± 0.100 (se)"
@test roundste(1049.0049,1.0) == "1049.00 ± 1.00 (se)"
@test roundste(1049.005,1.0) == "1049.01 ± 1.00 (se)"
@test roundste(0.,34567.) == "0 ± 34600 (se)"
@test roundste(1,NaN) == roundste(1.0,NaN) == "1 ± NaN (se)"
@test roundste(1.01,NaN) == "1.01 ± NaN (se)"
@test roundste(NaN,NaN) == roundste(NaN, 12.) == "NaN"

@test roundste(-100.,34567.) == "-100 ± 34600 (se)"
@test roundste(-0.,12345.) == "0 ± 12300 (se)"
@test roundste(-100.,4567.) == "-100 ± 4570 (se)"
@test roundste(-1099.,4567.) == "-1100 ± 4570 (se)"
@test roundste(-1049.,45649.) == "-1000 ± 45600 (se)"
@test roundste(-1049.,0.001) == "-1049.00000 ± 0.00100 (se)"
@test roundste(-1049.,0.1) == "-1049.000 ± 0.100 (se)"
@test roundste(-1049.0049,1.0) == "-1049.00 ± 1.00 (se)"
@test roundste(-1049.005,1.0) == "-1049.01 ± 1.00 (se)"
@test roundste(-0.,34567.) == "0 ± 34600 (se)"
@test roundste(-1,NaN) == roundste(-1.0,NaN) == "-1 ± NaN (se)"
@test roundste(-1.01,NaN) == "-1.01 ± NaN (se)"


# signif(1/10001,3) == 9.999999999999999e-5
#1//10001 == 0.000099990000999900009999000099990000...

@test roundste(1/10001,.001) == "0.00010 ± 0.00100 (se)"
@test roundste(1/10001,.0001) == "0.000100 ± 0.000100 (se)"
@test roundste(1/10001,.00001) == "0.0001000 ± 0.0000100 (se)"
@test roundste(1/10001,.000001) == "0.00009999 ± 0.00000100 (se)"
