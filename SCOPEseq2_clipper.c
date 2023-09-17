#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() 
{
	int i = 1, pos, clip = 0;
	char polya[] = "AAAAAAAA", polyn[] = "NNN", newline[500];
	char *line = NULL, *pa;
	size_t len = 0; 
	
	while (getline(&line, &len, stdin) != -1) {
		if (i == 1) {
			printf("%s",line);
			i = 2;
		}
		else if (i == 2) {
			pa = strstr(line, polya);
			if (pa == NULL) {
				printf("%s",line);
			}
			else {
				pos = pa - line;
				memset(newline,0,sizeof(newline));
				if (pos > 0) {
					strncpy(newline,line,pos);
					printf("%s\n",newline);
				}
				else {
					pos = 3;
					printf("%s\n",polyn);
				}
				clip++;
			}
			i = 3;
		}
		else if (i == 3) {
			printf("%s",line);
			i = 4;
		}
		else if (i == 4) {
			if (pa == NULL) {
				printf("%s",line);
				i = 1;
			}
			else {
				memset(newline,0,sizeof(newline));
				strncpy(newline,line,pos);
				printf("%s\n",newline);
			}
			i = 1;
		}
	}
	return clip;
}
