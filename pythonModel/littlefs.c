#include <stdio.h>
#include <stdint.h>
#include <string.h>

// Define before including lfs.h to disable malloc (critical for embedded)
#define LFS_NO_MALLOC

#include "lfs.h"

#define BLOCK_SIZE      4096
#define BLOCK_COUNT     4
#define CACHE_SIZE      256
#define LOOKAHEAD_SIZE  16
#define IMAGE_SIZE      (BLOCK_SIZE * BLOCK_COUNT)

uint8_t flash[IMAGE_SIZE];  // RAM-backed "flash"

// Static buffers for embedded (avoids malloc)
// Must be 4-byte aligned for ARM Cortex-M
static uint8_t read_buffer[CACHE_SIZE] __attribute__((aligned(4)));
static uint8_t prog_buffer[CACHE_SIZE] __attribute__((aligned(4)));
static uint8_t lookahead_buffer[LOOKAHEAD_SIZE] __attribute__((aligned(4)));
static uint8_t file_buffer[CACHE_SIZE] __attribute__((aligned(4)));  // For file operations
// -------------------- Block device functions --------------------
static int bd_read(const struct lfs_config *c, lfs_block_t block, lfs_off_t off,
                   void *buffer, lfs_size_t size) {
    memcpy(buffer, &flash[block * c->block_size + off], size);
    return 0;
}

static int bd_prog(const struct lfs_config *c, lfs_block_t block, lfs_off_t off,
                   const void *buffer, lfs_size_t size) {
    uint8_t *dst = &flash[block * c->block_size + off];
    const uint8_t *src = buffer;
    for (lfs_size_t i = 0; i < size; i++) dst[i] &= src[i]; // flash 1->0
    return 0;
}

static int bd_erase(const struct lfs_config *c, lfs_block_t block) {
    memset(&flash[block * c->block_size], 0xFF, c->block_size);
    return 0;
}

static int bd_sync(const struct lfs_config *c) { return 0; }

// -------------------- Main workflow --------------------
int main(void) {
    static struct lfs_config cfg = {0};  // Static to reduce stack usage
    static lfs_t lfs;                    // Static - lfs_t is large (~600 bytes)

    // Configure LittleFS
    cfg.read  = bd_read;
    cfg.prog  = bd_prog;
    cfg.erase = bd_erase;
    cfg.sync  = bd_sync;
    cfg.read_size      = 16;
    cfg.prog_size      = 16;
    cfg.block_size     = BLOCK_SIZE;
    cfg.block_count    = BLOCK_COUNT;
    cfg.block_cycles   = 500;  // Wear leveling cycles (required)
    cfg.cache_size     = CACHE_SIZE;
    cfg.lookahead_size = LOOKAHEAD_SIZE;

    // Provide static buffers (critical for embedded - avoids malloc)
    cfg.read_buffer      = read_buffer;
    cfg.prog_buffer      = prog_buffer;
    cfg.lookahead_buffer = lookahead_buffer;

    // Format the filesystem in RAM
    if (lfs_format(&lfs, &cfg) != 0) {
        printf("Failed to format LittleFS!\n");
        return 1;
    }
    printf("Filesystem formatted in RAM\n");

    // Mount it
    if (lfs_mount(&lfs, &cfg) != 0) {
        printf("Failed to mount LittleFS!\n");
        return 1;
    }
    printf("Filesystem mounted\n");

    // Create a CSV file
    static lfs_file_t file;          // Static - lfs_file_t is large
    static struct lfs_file_config file_cfg = {0};
    file_cfg.buffer = file_buffer;   // Static buffer for file ops
    
    if (lfs_file_opencfg(&lfs, &file, "csv1.csv", 
            LFS_O_WRONLY | LFS_O_CREAT, &file_cfg) != 0) {
        printf("Failed to create CSV file\n");
        return 1;
    }

    const char *csv_data = "header1,header2,header3\n1,2,3\n";
    lfs_file_write(&lfs, &file, csv_data, strlen(csv_data));
    lfs_file_close(&lfs, &file);
    printf("CSV file written\n");

    // Unmount filesystem
    lfs_unmount(&lfs);

    // Dump RAM image to a header
    FILE *f = fopen("lfs_image.h", "w");
    if (!f) return 1;

    fprintf(f, "unsigned char lfs_bin[%d] = {", IMAGE_SIZE);
    for (int i = 0; i < IMAGE_SIZE; i++) {
        if (i % 12 == 0) fprintf(f, "\n    ");
        fprintf(f, "0x%02X,", flash[i]);
    }
    fprintf(f, "\n};\n");
    fclose(f);

    printf("Header lfs_image.h written (%d bytes)\n", IMAGE_SIZE);
    return 0;
}
