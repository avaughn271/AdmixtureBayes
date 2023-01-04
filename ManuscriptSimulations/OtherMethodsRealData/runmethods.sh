treemix -i ArcticData.txt.gz -o output3 -root Yoruba -k 1000 -m 3 -global
treemix -i ArcticData.txt.gz -o output4 -root Yoruba -k 1000 -m 4 -global
treemix -i ArcticData.txt.gz -o output5 -root Yoruba -k 1000 -m 5 -global

orientagraph-i ArcticData.txt.gz -o orient3 -k 1000 -global -mlno -allmigs -root Yoruba -m 3
orientagraph -i ArcticData.txt.gz -o orient4 -k 1000 -global -mlno 1,2,3 -allmigs 1,2,3 -root Yoruba -m 4
orientagraph -i ArcticData.txt.gz -o orient5 -k 1000 -global  -mlno 1,2 -allmigs 1,2 -root Yoruba -m 5
