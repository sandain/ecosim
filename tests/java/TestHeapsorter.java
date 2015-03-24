
import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.Heapsorter;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class TestHeapsorter {

    @Before
    public void setup () {
        list = Arrays.asList (9, 8, 0, 2, 3, 1, 4, 5, 7, 6, -1);
        sorter = new Heapsorter<Integer> ();
    }

    @Test
    public void testReverseOrder () {
        sorter.sort (list);
        assertEquals ("Out of order.", 9, list.get (0).intValue ());
        assertEquals ("Out of order.", 8, list.get (1).intValue ());
        assertEquals ("Out of order.", 7, list.get (2).intValue ());
        assertEquals ("Out of order.", 6, list.get (3).intValue ());
        assertEquals ("Out of order.", 5, list.get (4).intValue ());
        assertEquals ("Out of order.", 4, list.get (5).intValue ());
        assertEquals ("Out of order.", 3, list.get (6).intValue ());
        assertEquals ("Out of order.", 2, list.get (7).intValue ());
        assertEquals ("Out of order.", 1, list.get (8).intValue ());
        assertEquals ("Out of order.", 0, list.get (9).intValue ());
        assertEquals ("Out of order.", -1, list.get (10).intValue ());
    }

    @Test
    public void testForwardOrder () {
        sorter.sort (list, new Comparator<Integer> () {
            @Override
            public int compare (Integer a, Integer b) {
                return b.compareTo (a);
            }
        });
        assertEquals ("Out of order.", -1, list.get (0).intValue ());
        assertEquals ("Out of order.", 0, list.get (1).intValue ());
        assertEquals ("Out of order.", 1, list.get (2).intValue ());
        assertEquals ("Out of order.", 2, list.get (3).intValue ());
        assertEquals ("Out of order.", 3, list.get (4).intValue ());
        assertEquals ("Out of order.", 4, list.get (5).intValue ());
        assertEquals ("Out of order.", 5, list.get (6).intValue ());
        assertEquals ("Out of order.", 6, list.get (7).intValue ());
        assertEquals ("Out of order.", 7, list.get (8).intValue ());
        assertEquals ("Out of order.", 8, list.get (9).intValue ());
        assertEquals ("Out of order.", 9, list.get (10).intValue ());
    }

    private List<Integer> list;
    private Heapsorter<Integer> sorter;

}
