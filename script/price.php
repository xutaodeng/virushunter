<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
 <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
  <title>Price</title>

  <script src="./tree/jquery/jquery.js" type="text/javascript"></script>
  <script src="./tree/jquery/jquery-ui.custom.js" type="text/javascript"></script>
  <script src="./tree/jquery/jquery.cookie.js" type="text/javascript"></script>

  <link href="./tree/src/skin/ui.dynatree.css" rel="stylesheet" type="text/css"> 
  <script src="./tree/src/jquery.dynatree.js" type="text/javascript"></script>
<?php
 function startsWith($haystack, $needle)
{
    return !strncmp($haystack, $needle, strlen($needle));
}

function endsWith($haystack, $needle)
{
    $length = strlen($needle);
    if ($length == 0) {
        return true;
    }
    return (substr($haystack, -$length) === $needle);
}
?>
 <script type="text/javascript">
  db='all'
  $(function(){
      $("#tree").dynatree({
      onActivate: function(node) {
	    var db = node.data.title;
	    document.forms["myForm"]["db"].value=db;
        // alert("You activated " + document.forms["myForm"]["db"].value);
        
      },
      children: [
	{title: "All", isFolder: true, key: "All",
	  children:[
<?php
	 $handle=opendir(".");
	 $allout=array();
	 
	while ($file = readdir($handle)) 
	{
		if (is_dir($file) && $file!='' && $file!='.' && $file!='..' && $file!='DataTables-1.9.4' && $file!='tree') 
		{
			$children=array();
			$parent="{title:\"NULL\", isFolder: true, children:[";;
			if (is_dir($file.'/blast'))
			{
				if ($handle2=opendir($file.'/blast'))
				{	$c=0;
					while ($file2 = readdir($handle2)) 
					{
						if (endsWith($file2, ".nhr") or endsWith($file2, ".nal")) 
						{
							$fid = substr($file2, 0, strpos($file2, "_blastdb"));
							// if (endsWith($fid, "_1") or endsWith($fid, "_2")){
								// continue;
							// }
							//{title: "Folder 2", isFolder: true, children: [{title: "Sub-item 2.1"},{title: "Sub-item 2.2"}]},
							if ($fid==$file){
								//echo 'parent? '.$file.'\n';
								$parent="{title:\"".$fid."\", isFolder: true, children:[";
							}
							else {
								$child='{title: '."\"".$file.'$'.$fid."\"}";
								if (!in_array($child, $children)){
									array_push($children,  $child);
									$c=$c+1;
								}
							}
						}
					}
					if ($c>0){
						$output=$parent.implode(",", $children)."]}\n";
						array_push($allout, $output);
					}
				}
			}
		}
	}
	echo implode(",", $allout);
 ?>
 	  ]
	}
      ]
    });
  });
 </script>

<script type="text/javascript">
function clearForms()
{
    document.getElementById('pleaseWait').reset();
}
function validateForm()
{

var x=document.forms["myForm"]["query"].value;
	if (x==null || x=="")
	  {
	  alert("You have empty query sequence!");
	  return false;
	  }
var x=document.forms["myForm"]["db"].value;
	if (x.indexOf('$') == -1){
	alert("select your database at barcode level!");
	return false;
	}
	document.getElementById('pleaseWait').style.display='block';
	return true;
}
</script>
</head>

  <DIV ID="pleaseWait" STYLE="position:absolute;z-index:5;top:30%;left:35%; display: none">
	<IMG SRC="wait.gif" WIDTH=600 HEIGHT=300>
		  <!-- <IMG SRC="wait.gif" BORDER=1 WIDTH=75 HEIGHT=15><BR><BR>  -->
  </DIV> 
 <body onLoad="clearForms()" onUnload="clearForms()">
  
<form name="myForm" action="price_run.php" method="post" onsubmit="return validateForm()">

Paste Your Seed Sequence: <br>
 <TEXTAREA name="query" rows="10" cols="80"></TEXTAREA><br>
<br>
Select Run to Price Against: <br>
 <div id="tree"> </div>

<br>
<INPUT type="hidden" name="db" value="All">
<FONT COLOR="red">Warning: Price consumes a lot of system resources. Always run ONE job at a time and it will take hours!!!</FONT>
<br>
<input type="submit" name = "submit" value="Run Price"  style="background-color:lightgreen"> <br>

</form> 
</body>
</html>