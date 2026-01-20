#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

int main(void) {
    int fp = open("ytrue_UKF_5hr_10s.csv", O_RDONLY);
    if (fp == -1) {
        perror("open");
        return 1;
    }

    char buf[256];
    char line[512];
    int line_len = 0;
    int n;

    while ((n = read(fp, buf, sizeof(buf))) > 0) {
        for (int i = 0; i < n; i++) {
            if (buf[i] == '\n') {
                line[line_len] = '\0';
                printf("%s\n", line);
                line_len = 0;
            } else {
                line[line_len++] = buf[i];
            }
        }
    }

    close(fp);
    return 0;
}
