Image Alignment Tool
(c) GreatAttractor (WielkiAtraktor)
v. 0.5
2014/05/22

Dozwolone rozpowszechnianie i wykorzystanie w dowolny spos�b bez ogranicze�.


1. Wprowadzenie
  1.1 Sk�adnia w linii polece�
  1.2 Pobieranie
  1.3 Uwagi
    1.3.1 Konwersja plik�w
    1.3.2 Uruchamianie bez u�ycia wiersza polece� (konsoli)
2. Budowanie ze �r�de�
2.1 Budowanie za pomoc� MinGW
2.2 Budowanie za pomoc� Microsoft C++ (z Windows SDK b�d� Visual Studio/Visual C++)
2.3 Budowanie za pomoc� GCC (Linux i pokrewne, MinGW+MSYS)
3 Historia zmian

----------------------------------------
1. Wprowadzenie

Imgalt jest narz�dziem do wyr�wnywania sekwencji obraz�w, kt�re mog� by� znacznie przesuni�te wzgl�dem siebie w przypadkowy spos�b. Do obr�bki przyjmowane s� pliki BMP (8- lub 24-bitowe) i TIFF (8 lub 16 bit�w na kana�, RGB lub w odcieniach szaro�ci), rozmiary mog� si� r�ni�. Po wyr�wnaniu pliki wynikowe zapisywane s� w tym samym formacie, o jednakowej wielko�ci r�wnej wielko�ci prostok�ta b�d�cego przeci�ciem (najwi�kszym wsp�lnym obszarem) wszystkich obraz�w wej�ciowych po ich wyr�wnaniu. (Mo�na te� zapisa� obrazy powi�kszone do obejmuj�cego je wszystkie prostok�ta u�ywaj�c parametru --no-crop).

Z powodu interpolacji niezb�dnej do wyr�wnywania subpikselowego, 8-bitowe wej�ciowe pliki BMP (z palet�) s� zapisywane jako 24-bitowe BMP RGB. (Nie ma to miejsca, gdy podany jest parametr "--no-subpixel").

Przyk�adowy scenariusz u�ycia to przygotowanie klatek animacji w osobnym katalogu (folderze), przekazanie ich do imgalt, otwarcie wyr�wnanych klatek w GIMPie (File->Open as Layers...), obr�bka i przyci�cie (zaznaczenie obszaru, nast�pnie Image->Crop to Selection) i wreszcie eksport jako animowany GIF (File->Export, wyb�r formatu GIF, zaznaczenie "As Animation").

  ----------------------------------------
  1.1 Sk�adnia w linii polece�
  
Og�lna posta� wywo�ania:
  
  imgalt [--verbose] [--threads <liczba>] [--no-crop] [--no-subpixel] [--output-dir <output directory>] [[--input-dir] <input directory>]
  
Katalog wej�ciowy jest parametrem domy�lnym, wi�c mo�na pomin�� "--input-dir". Przyk�ady:

  imgalt c:\astro\solar
  
(wyr�wnuje wszystkie pliki w c:\astro\solar i zapisuje wyr�wnane w katalogu bie��cym).

  imgalt
  
(wyr�wnuje wszystkie pliki w katalogu bie��cym i zapisuje wyr�wnane w tym samym miejscu).

Imgalt wyr�wnuje wszystkie pliki BMP i TIFF (posortowane po nazwie) w katalogu wej�ciowym (domy�lnie: bie��cy) i zapisuje je z przyrostkiem "_aligned" w katalogu wyj�ciowym (domy�lnie: bie��cy).

Parametr "--verbose" w��cza wypisywanie dodatkowych informacji w trakcie wyr�wnywania.

Parametr "--threads" ustawia liczb� w�tk�w roboczych (domy�lnie: liczba wykrytych procesor�w logicznych).

Parametr "--no-crop" powoduje powi�kszenie obraz�w wyj�ciowych do rozmiaru wsp�lnego prostok�ta otaczaj�cego zamiast obcinania ich do najwi�kszego wsp�lnego obszaru.

Parametr "--no-subpixel" wy��cza wyr�wnywanie subpikselowe. Zastosowanie go mo�e spowodowa� bardziej widoczny "dryf" obrazu wyj�ciowego lub nieznaczne 1-pikselowe skoki. Ponadto zapis obraz�w wyj�ciowych b�dzie szybszy (ale nie samo wyr�wnywanie).

Je�li katalog wej�ciowy lub docelowy zawieraj� spacje, nale�y otoczy� �cie�ki znakami cudzys�owu, np.:

  imgalt "c:\astro\Czerwiec 12" --output-dir "c:\astro\Czerwiec 12\aligned"

  ----------------------------------------
  1.2 Pobieranie
  
Najnowsz� wersj� mo�na pobra� z:

http://stargazerslounge.com/blog/1400/entry-1654-imgalt/
http://astropolis.pl/topic/44806-narzedzie-do-automatycznego-wyrownywania-klatek-animacji-slonecznych/#entry534928
  
  ----------------------------------------
  1.3 Uwagi

    ----------------------------------------
    1.3.1 Konwersja plik�w
    
Pliki mo�na szybko skonwertowa� zbiorowo np. w IrfanView (File->Batch Conversion/Rename...). Animowany GIF mo�na rozbi� na klatki np. w programie VirtualDub (otworzy� plik GIF, potem File->Export->Image sequence...).
  
    ----------------------------------------
    1.3.2 Uruchamianie bez u�ycia wiersza polece� (konsoli)
    
Jako �e uruchomienie imgalt bez parametr�w oznacza przyj�cie katalogu (folderu) bie��cego jako wej�ciowego, mo�na zrobi�, co nast�puje:
  - skopiowa� pliki przeznaczone do wyr�wnania do katalogu, gdzie znajduje si� imgalt.exe
  - uruchomi� imgalt.exe
Wyr�wnanie pliki wyj�ciowe pojawi� si� w tym samym katalogu.

Kolejn� metod� (pod MS Windows, cho� mo�e zadzia�a� r�wnie� w innych �rodowiskach graficznych) jest przeci�gni�cie katalogu z plikami wej�ciowymi i upuszczenie go na ikon� imgalt.exe. Wyr�wnane pliki pojawi� si� w katalogu, kt�ry by� w tym czasie bie��cy (np. w Windows, je�li przeci�gni�ty zosta� "d:\astro\sekwencja1", pliki wynikowe powinny pojawi� si� w "d:\astro").
  
----------------------------------------
2. Budowanie ze �r�de�

Kod �r�d�owy imgalt mo�na swobodnie rozpowszechnia� i wykorzystywa� w dowolnym celu. Do wybudowania ze �r�de� potrzebne s� biblioteki Boost w wersji 1.53.0 lub nowszej (starsze mog� dzia�a� po niewielkich zmianach w kodzie imgalt). Wielow�tkowo�� wymaga kompilatora obs�uguj�cego OpenMP (np. GCC 4.2 lub nowsze, MS Visual C++ 2008 lub nowsze (wersje p�atne), MS Visual C++ 2012 Express lub nowsze).

  ----------------------------------------
  2.1 Budowanie za pomoc� MinGW

Otworzy� Makefile.mingw i ustawi� zmienne BOOST_INCLUDE i BOOST_LIBS na odpowiednie dla u�ywanej instalacji Boost. Upewni� si�, �e narz�dzia MinGW s� na �cie�ce (np. wywo�a� "set PATH=%PATH%;c:\mingw\bin") i wykona�:

    mingw32-make -f Makefile.mingw
    
  ----------------------------------------    
  2.2 Budowanie za pomoc� Microsoft C++ (z Windows SDK b�d� Visual Studio/Visual C++)
  
Otworzy� Makefile.msvc i ustawi� zmienne BOOST_INCLUDE i BOOST_LIBS na odpowiednie dla u�ywanej instalacji Boost. Przeczyta� notatk� o obs�udze OpenMP. Upewni� si�, �e narz�dzia MS s� na �cie�ce (np. wywo�a� "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\vcvars32.bat") i wykona�:

    nmake -f Makefile.msvc
    
  ----------------------------------------        
  2.3 Budowanie za pomoc� GCC (Linux i pokrewne, MinGW+MSYS)
  
Upewni� si�, �e binarne i nag��wkowe biblioteki Boost s� zainstalowane (np. pod Linuksem zainstalowa� pakiety "boost", "boost-devel") i wykona�:

    make -f Makefile.gcc
    
----------------------------------------            
3. Historia zmian

0.5 (2014/05/22)
    Nowe funkcje:
      - wyr�wnywanie subpikselowe

0.4.1 (2014/05/05)
    Nowe funkcje:
      - obs�uga plik�w TIFF

0.4 (2014/05/02)
    Nowe funkcje:
      - domy�lne obcinanie do najwi�kszego wsp�lnego obszaru
    Poprawki b��d�w:
      - poprawiona translacja przy zapisie wyr�wnanych obraz�w, powodowa�a chaotyczne 1-pikselowe skoki

0.3 (2014/04/25)
    Nowe funkcje:
      - przyjmowanie r�wnie� kolorowych obraz�w wej�ciowych
      - zapis obraz�w wyj�ciowych w tym samym formacie kolor�w co wej�ciowe
    Poprawki b��d�w:
      - poprawiony zapis BMP o d�ugo�ci linii nie b�d�cej wielokrotno�ci� 4