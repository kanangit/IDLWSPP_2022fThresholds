;converts units from pixels to mm, shift the origin to zero if needed
FUNCTION px_to_mm, pos_px, shift_px, scale_mm_in_px

pos_mm = (DOUBLE(pos_px) - DOUBLE(shift_px)) * DOUBLE(scale_mm_in_px)

RETURN, pos_mm
END