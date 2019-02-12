"c:\Program Files (x86)\Resource Hacker\ResHacker.exe" -addoverwrite dke.exe, dkenew.exe, dke.ico, ICONGROUP,MAINICON,0

move dkenew.exe dke.exe

"c:\Program Files (x86)\Resource Hacker\ResHacker.exe" -addoverwrite dke_preprocess_dicom.exe, dke_preprocess_dicomnew.exe, dke.ico, ICONGROUP,MAINICON,0

move dke_preprocess_dicomnew.exe dke_preprocess_dicom.exe

"c:\Program Files (x86)\Resource Hacker\ResHacker.exe" -addoverwrite map_interpolate.exe, map_interpolatenew.exe, dke.ico, ICONGROUP,MAINICON,0

move map_interpolatenew.exe map_interpolate.exe
