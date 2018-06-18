# initialize the tfTarget package

.onAttach<- function(libname, pkgName)
{
	assignInNamespace("background.check", background.check,ns="rtfbsdb")
	assignInNamespace("tfbs_enrichmentTest", tfbs_enrichmentTest,ns="rtfbsdb")
}