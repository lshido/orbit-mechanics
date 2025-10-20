from great_tables import GT, md, html
from great_tables.data import islands

islands_mini = islands.head(10)
# Create a display table showing ten of the largest islands in the world
gt_tbl = GT(islands_mini)

# Show the output table
gt_tbl.show()