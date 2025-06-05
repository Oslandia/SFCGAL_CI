// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2025 Oslandia
// Copyright (c) 2024-2025, SFCGAL team.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SFCGAL/capi/sfcgal_c.h"

/* This is a very simple c program to check that the C API compiles.
 */
int
main(void)
{
  sfcgal_init();
  char               wkt_str[] = "LINESTRING (0.0 0.0,2.0 0.0,1.0 1.0)";
  sfcgal_geometry_t *geometry  = sfcgal_io_read_wkt(wkt_str, strlen(wkt_str));

  char  *wkt_buffer;
  size_t wkt_len;
  sfcgal_geometry_as_text_decim(geometry, 1, &wkt_buffer, &wkt_len);

  if (strcmp(wkt_buffer, wkt_str)) {
    fprintf(stderr, "WKT do not match\n");
    sfcgal_geometry_delete(geometry);
    sfcgal_free_buffer(wkt_buffer);
    return EXIT_FAILURE;
  }

  printf("WKT match\n");
  sfcgal_geometry_delete(geometry);
  sfcgal_free_buffer(wkt_buffer);

  return EXIT_SUCCESS;
}
