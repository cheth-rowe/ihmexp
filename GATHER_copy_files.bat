REM 2019 GATHER batch process
j:
cd \Project\GBD\GBD_2019_Code
c:
cd \data\repos\IT\ihmexp

REM --remove all existing C:repo files first---
del /s /q \data\repos\IT\ihmexp\gbd_2019\appendices\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\cod_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\fert_pop_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\mortality_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\nonfatal_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\risk_factors_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\sdgs_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\shared_code\*
del /s /q \data\repos\IT\ihmexp\gbd_2019\uhc_code\*

xcopy j:appendices              C:\data\repos\IT\ihmexp\gbd_2019\appendices /s /y
xcopy j:cod_code                C:\data\repos\IT\ihmexp\gbd_2019\cod_code /s /y
xcopy j:fert_pop_code           C:\data\repos\IT\ihmexp\gbd_2019\fert_pop_code /s /y
xcopy j:mortality_code          C:\data\repos\IT\ihmexp\gbd_2019\mortality_code /s /y
xcopy j:nonfatal_code           C:\data\repos\IT\ihmexp\gbd_2019\nonfatal_code /s /y
xcopy j:risk_factors_code       C:\data\repos\IT\ihmexp\gbd_2019\risk_factors_code /s /y
xcopy j:sdgs_code               C:\data\repos\IT\ihmexp\gbd_2019\sdgs_code /s /y
xcopy j:shared_code             C:\data\repos\IT\ihmexp\gbd_2019\shared_code /s /y
xcopy j:uhc_code                C:\data\repos\IT\ihmexp\gbd_2019\uhc_code /s /y
