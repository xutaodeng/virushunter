<html>
<body>

Your seeds is: <?php echo $_POST['query']; ?><br />
Your barcode is: <?php echo $_POST['db']; ?><br />
<?php 
set_time_limit(0);
$stamp = date("Y-m-d-H-i-s");
//$stamp = '2013-02-05-05-29-33';
$seeds = 'tmp\\'.$stamp."_seeds.fa";
$fh = fopen($seeds, 'w') or die("can't open file");
$stringData = trim($_POST['query']);
$price="C:\\cygwin\\bin\\PriceTI.exe";
$db=$_POST['db'];
//echo $stringData;
if ($stringData[0]!='>'){fwrite($fh, "> query\n");}
fwrite($fh, $stringData);
fclose($fh);
?>
<br />

<?php
//$wd=getcwd();
$wd='/cygdrive/c/wamp/www';
echo '<pre>';
$path=explode('$', $db);
$library=$path[0];
$barcode=$path[1];
echo $library;
echo '<br>'.$barcode.'<br>';
$dbfile1=$wd.'/'.$library.'/price/'.$barcode.'_1.fq';
$dbfile2=$wd.'/'.$library.'/price/'.$barcode.'_2.fq';
$outfa=$wd.'/tmp/'.$library.$stamp.'.fa';

$cmd2 =$price.' -icf '.$wd.'/'.$seeds.' 1 1 5'.' -fp  '.$dbfile1.' '.$dbfile2.' 300 -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -o '.$outfa;

// //PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  
#$cmd2='C:\\cygwin\\bin\\PriceTI.exe > error.log 2>&1';
//$cmd2 = 'C:\\bin\\PriceTI.exe > error.log 2>&1';
//$cmd2='notepad.exe';
echo '<br>'.$cmd2.'<br>';
$retval='';
$last_line = system($cmd2, $retval);
//$last_line = shell_exec($cmd2);
// Printing additional info
 echo '
 </pre>
  <hr />Last line of the output: ' . $last_line . '
 <hr />Return value: ' . $retval.'
 <br />';

echo '<br>';
for ($i=1; $i<=30; $i++)
{
	//$filename ='./'.$library.'/price/'.$stamp.'.cycle'.$i.'.fa';
	$filename ='./tmp/'.$library.$stamp.'.cycle'.$i.'.fa';
	//echo $filename.'<br>';
	if (file_exists($filename)){
		echo '<a href="'.$filename.'"</a>'.'output_cycle'.$i.'<br>'; 
	}
	// else{
		// break;
	// }
}

// $file =  "price30.fa"; //Path to your *.txt file
// $contents = file($file);
// $string = implode($contents);
// echo $string;
?>

</body>
</html>