package Homyaki::Geo_Maps::Conversion::Pulkovo;

use Math::Complex;

use constant PI => 3.14159265358979; # Число Пи
use constant RO => 206264.8062;      # Число угловых секунд в радиане

# Эллипсоид Красовского
use constant AP  => 6378245;            # Большая полуось
use constant AlP => 1 / 298.3;          # Сжатие
use constant E2P => 2 * AlP - AlP ** 2; # Квадрат эксцентриситета
                                                                                                                     
# Эллипсоид WGS84 (GRS80, эти два эллипсоида сходны по большинству параметров)
use constant AW  => 6378137;            # Большая полуось
use constant AlW => 1 / 298.257223563;  # Сжатие
use constant E2W => 2 * AlW - AlW ** 2; # Квадрат эксцентриситета

# Вспомогательные значения для преобразования эллипсоидов
use constant A  => (AP + AW) / 2;
use constant E2  => (E2P + E2W) / 2;
use constant DA  => AW - AP;
use constant DE2 => E2W - E2P;

# Линейные элементы трансформирования, в метрах
use constant DX => 23.92;
use constant DY => -141.27;
use constant DZ => -80.9;

# Угловые элементы трансформирования, в секундах
use constant WX => 0;
use constant WY => 0;
use constant WZ => 0;

# Дифференциальное различие масштабов
use constant MS => 0;

sub convert {
	my $class = shift;
	my %h = @_;

	my $lat = $h{lat};
	my $lng = $h{lng};

	return {
		lat => SK42_WGS84_Lat($lat, $lng, 0),
		lng => SK42_WGS84_Long($lat, $lng, 0),
	};
}

sub WGS84_SK42_Lat {
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	return $bd - db($bd, $ld, $h) / 3600;
}

sub SK42_WGS84_Lat {
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	return $bd + db($bd, $ld, $h) / 3600;
}

sub WGS84_SK42_Long {                                                                                                                                                                                        
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	return $ld - dl($bd, $ld, $h) / 3600;
}

sub SK42_WGS84_Long {
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	return $ld + dl($bd, $ld, $h) / 3600;
}

sub db{
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	$b = $bd * &PI / 180;
	$l = $ld * &PI / 180;
	$m = &A * (1 - &E2) / (1 - &E2 * sin($b) ** 2) ** 1.5;
	$n = &A * (1 - &E2 * sin($b) ** 2) ** -0.5;

	return &RO / ($m + $h) * ($n / &A * &E2 * sin($b) * cos($b) * &DA 
	 + ($n ** 2 / &A ** 2 + 1) * $n * sin($b) * cos($b) * &DE2 / 2 
	 - (&DX * cos($l) + &DY * sin($l)) * sin($b) + &DZ * cos($b)) 
	 - &WX * sin($l) * (1 + &E2 * cos(2 * $b)) 
	 + &WY * cos($l) * (1 + &E2 * cos(2 * $b)) 
	 - &RO * &MS * &E2 * sin($b) * cos($b);

}

sub dl{
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	$b = $bd * &PI / 180;
	$l = $ld * &PI / 180;
	$n = &A * (1 - &E2 * sin($b) ** 2) ** -0.5;

	return &RO / (($n + $h) * cos($b)) * (-&DX * sin($l) + &DY * cos($l)) 
	 + tan($b) * (1 - &E2) * (&WX * cos($l) + &WY * sin($l)) - &WZ;
}

sub WGS84Alt{
	my $bd = shift;
	my $ld = shift;
	my $h  = shift;

	my $b = $bd * &PI / 180;
	my $l = $ld * &PI / 180;
	my $n = &A * (1 - &E2 * sin($b) ** 2) ** -0.5;

	my $dh = -&A / $n * &DA + $n * sin($b) ** 2 * &DE2 / 2
	 + (&DX * cos($l) + &DY * sin($l)) * cos($b) + &DZ * sin($b)
	 - $n * &E2 * sin($b) * cos($b) * (&WX / &RO * sin($l) - &WY / &RO * cos($l)) 
	 + (&A ** 2 / $n + $h) * &MS;

	return $h + $dh;
}

1;