#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>

#define READ_BUF 512
#define LINE_BUF 2048
#define MAX_COLS 16   // SoC + 15 temperatures

int main(void)
{
    int fd = open("Sample_OCV_SoC_1.0pct_10deg.csv", O_RDONLY);
    if (fd < 0) {
        perror("open");
        return 1;
    }

    char rbuf[READ_BUF];
    char line[LINE_BUF];
    int line_len = 0;
    int header_skipped = 0;

    ssize_t n;

    while ((n = read(fd, rbuf, sizeof(rbuf))) > 0) {
        for (ssize_t i = 0; i < n; i++) {

            if (rbuf[i] == '\n') {
                line[line_len] = '\0';

                /* Skip header line */
                if (!header_skipped) {
                    header_skipped = 1;
                    line_len = 0;
                    continue;
                }

                /* ---- Parse one CSV row ---- */
                double values[MAX_COLS];
                int col = 0;

                char *p = line;
                char *end;

                while (*p && col < MAX_COLS) {
                    values[col++] = strtod(p, &end);
                    if (p == end) break;
                    p = end;
                    if (*p == ',') p++;
                }

                /* Example usage */
                printf("SoC=%0.f  OCV@20C=%f  OCV@40C=%f\n",
                       values[0], values[7], values[9]);

                line_len = 0;
            }
            else {
                if (line_len < LINE_BUF - 1) {
                    line[line_len++] = rbuf[i];
                }
            }
        }
    }

    close(fd);
    return 0;
}
