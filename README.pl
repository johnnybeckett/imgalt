Image Alignment Tool
(c) GreatAttractor (WielkiAtraktor)
v. 0.5
2014/05/22

Dozwolone rozpowszechnianie i wykorzystanie w dowolny sposób bez ograniczeñ.


1. Wprowadzenie
  1.1 Sk³adnia w linii poleceñ
  1.2 Pobieranie
  1.3 Uwagi
    1.3.1 Konwersja plików
    1.3.2 Uruchamianie bez u¿ycia wiersza poleceñ (konsoli)
2. Budowanie ze Ÿróde³
2.1 Budowanie za pomoc¹ MinGW
2.2 Budowanie za pomoc¹ Microsoft C++ (z Windows SDK b¹dŸ Visual Studio/Visual C++)
2.3 Budowanie za pomoc¹ GCC (Linux i pokrewne, MinGW+MSYS)
3 Historia zmian

----------------------------------------
1. Wprowadzenie

Imgalt jest narzêdziem do wyrównywania sekwencji obrazów, które mog¹ byæ znacznie przesuniête wzglêdem siebie w przypadkowy sposób. Do obróbki przyjmowane s¹ pliki BMP (8- lub 24-bitowe) i TIFF (8 lub 16 bitów na kana³, RGB lub w odcieniach szaroœci), rozmiary mog¹ siê ró¿niæ. Po wyrównaniu pliki wynikowe zapisywane s¹ w tym samym formacie, o jednakowej wielkoœci równej wielkoœci prostok¹ta bêd¹cego przeciêciem (najwiêkszym wspólnym obszarem) wszystkich obrazów wejœciowych po ich wyrównaniu. (Mo¿na te¿ zapisaæ obrazy powiêkszone do obejmuj¹cego je wszystkie prostok¹ta u¿ywaj¹c parametru --no-crop).

Z powodu interpolacji niezbêdnej do wyrównywania subpikselowego, 8-bitowe wejœciowe pliki BMP (z palet¹) s¹ zapisywane jako 24-bitowe BMP RGB. (Nie ma to miejsca, gdy podany jest parametr "--no-subpixel").

Przyk³adowy scenariusz u¿ycia to przygotowanie klatek animacji w osobnym katalogu (folderze), przekazanie ich do imgalt, otwarcie wyrównanych klatek w GIMPie (File->Open as Layers...), obróbka i przyciêcie (zaznaczenie obszaru, nastêpnie Image->Crop to Selection) i wreszcie eksport jako animowany GIF (File->Export, wybór formatu GIF, zaznaczenie "As Animation").

  ----------------------------------------
  1.1 Sk³adnia w linii poleceñ
  
Ogólna postaæ wywo³ania:
  
  imgalt [--verbose] [--threads <liczba>] [--no-crop] [--no-subpixel] [--output-dir <output directory>] [[--input-dir] <input directory>]
  
Katalog wejœciowy jest parametrem domyœlnym, wiêc mo¿na pomin¹æ "--input-dir". Przyk³ady:

  imgalt c:\astro\solar
  
(wyrównuje wszystkie pliki w c:\astro\solar i zapisuje wyrównane w katalogu bie¿¹cym).

  imgalt
  
(wyrównuje wszystkie pliki w katalogu bie¿¹cym i zapisuje wyrównane w tym samym miejscu).

Imgalt wyrównuje wszystkie pliki BMP i TIFF (posortowane po nazwie) w katalogu wejœciowym (domyœlnie: bie¿¹cy) i zapisuje je z przyrostkiem "_aligned" w katalogu wyjœciowym (domyœlnie: bie¿¹cy).

Parametr "--verbose" w³¹cza wypisywanie dodatkowych informacji w trakcie wyrównywania.

Parametr "--threads" ustawia liczbê w¹tków roboczych (domyœlnie: liczba wykrytych procesorów logicznych).

Parametr "--no-crop" powoduje powiêkszenie obrazów wyjœciowych do rozmiaru wspólnego prostok¹ta otaczaj¹cego zamiast obcinania ich do najwiêkszego wspólnego obszaru.

Parametr "--no-subpixel" wy³¹cza wyrównywanie subpikselowe. Zastosowanie go mo¿e spowodowaæ bardziej widoczny "dryf" obrazu wyjœciowego lub nieznaczne 1-pikselowe skoki. Ponadto zapis obrazów wyjœciowych bêdzie szybszy (ale nie samo wyrównywanie).

Jeœli katalog wejœciowy lub docelowy zawieraj¹ spacje, nale¿y otoczyæ œcie¿ki znakami cudzys³owu, np.:

  imgalt "c:\astro\Czerwiec 12" --output-dir "c:\astro\Czerwiec 12\aligned"

  ----------------------------------------
  1.2 Pobieranie
  
Najnowsz¹ wersjê mo¿na pobraæ z:

http://stargazerslounge.com/blog/1400/entry-1654-imgalt/
http://astropolis.pl/topic/44806-narzedzie-do-automatycznego-wyrownywania-klatek-animacji-slonecznych/#entry534928
  
  ----------------------------------------
  1.3 Uwagi

    ----------------------------------------
    1.3.1 Konwersja plików
    
Pliki mo¿na szybko skonwertowaæ zbiorowo np. w IrfanView (File->Batch Conversion/Rename...). Animowany GIF mo¿na rozbiæ na klatki np. w programie VirtualDub (otworzyæ plik GIF, potem File->Export->Image sequence...).
  
    ----------------------------------------
    1.3.2 Uruchamianie bez u¿ycia wiersza poleceñ (konsoli)
    
Jako ¿e uruchomienie imgalt bez parametrów oznacza przyjêcie katalogu (folderu) bie¿¹cego jako wejœciowego, mo¿na zrobiæ, co nastêpuje:
  - skopiowaæ pliki przeznaczone do wyrównania do katalogu, gdzie znajduje siê imgalt.exe
  - uruchomiæ imgalt.exe
Wyrównanie pliki wyjœciowe pojawi¹ siê w tym samym katalogu.

Kolejn¹ metod¹ (pod MS Windows, choæ mo¿e zadzia³aæ równie¿ w innych œrodowiskach graficznych) jest przeci¹gniêcie katalogu z plikami wejœciowymi i upuszczenie go na ikonê imgalt.exe. Wyrównane pliki pojawi¹ siê w katalogu, który by³ w tym czasie bie¿¹cy (np. w Windows, jeœli przeci¹gniêty zosta³ "d:\astro\sekwencja1", pliki wynikowe powinny pojawiæ siê w "d:\astro").
  
----------------------------------------
2. Budowanie ze Ÿróde³

Kod Ÿród³owy imgalt mo¿na swobodnie rozpowszechniaæ i wykorzystywaæ w dowolnym celu. Do wybudowania ze Ÿróde³ potrzebne s¹ biblioteki Boost w wersji 1.53.0 lub nowszej (starsze mog¹ dzia³aæ po niewielkich zmianach w kodzie imgalt). Wielow¹tkowoœæ wymaga kompilatora obs³uguj¹cego OpenMP (np. GCC 4.2 lub nowsze, MS Visual C++ 2008 lub nowsze (wersje p³atne), MS Visual C++ 2012 Express lub nowsze).

  ----------------------------------------
  2.1 Budowanie za pomoc¹ MinGW

Otworzyæ Makefile.mingw i ustawiæ zmienne BOOST_INCLUDE i BOOST_LIBS na odpowiednie dla u¿ywanej instalacji Boost. Upewniæ siê, ¿e narzêdzia MinGW s¹ na œcie¿ce (np. wywo³aæ "set PATH=%PATH%;c:\mingw\bin") i wykonaæ:

    mingw32-make -f Makefile.mingw
    
  ----------------------------------------    
  2.2 Budowanie za pomoc¹ Microsoft C++ (z Windows SDK b¹dŸ Visual Studio/Visual C++)
  
Otworzyæ Makefile.msvc i ustawiæ zmienne BOOST_INCLUDE i BOOST_LIBS na odpowiednie dla u¿ywanej instalacji Boost. Przeczytaæ notatkê o obs³udze OpenMP. Upewniæ siê, ¿e narzêdzia MS s¹ na œcie¿ce (np. wywo³aæ "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\vcvars32.bat") i wykonaæ:

    nmake -f Makefile.msvc
    
  ----------------------------------------        
  2.3 Budowanie za pomoc¹ GCC (Linux i pokrewne, MinGW+MSYS)
  
Upewniæ siê, ¿e binarne i nag³ówkowe biblioteki Boost s¹ zainstalowane (np. pod Linuksem zainstalowaæ pakiety "boost", "boost-devel") i wykonaæ:

    make -f Makefile.gcc
    
----------------------------------------            
3. Historia zmian

0.5 (2014/05/22)
    Nowe funkcje:
      - wyrównywanie subpikselowe

0.4.1 (2014/05/05)
    Nowe funkcje:
      - obs³uga plików TIFF

0.4 (2014/05/02)
    Nowe funkcje:
      - domyœlne obcinanie do najwiêkszego wspólnego obszaru
    Poprawki b³êdów:
      - poprawiona translacja przy zapisie wyrównanych obrazów, powodowa³a chaotyczne 1-pikselowe skoki

0.3 (2014/04/25)
    Nowe funkcje:
      - przyjmowanie równie¿ kolorowych obrazów wejœciowych
      - zapis obrazów wyjœciowych w tym samym formacie kolorów co wejœciowe
    Poprawki b³êdów:
      - poprawiony zapis BMP o d³ugoœci linii nie bêd¹cej wielokrotnoœci¹ 4